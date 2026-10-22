# =============================================================================
# gasampling: sampling_frame()
#
# Takes SSU-level census data and produces:
#   (1) a strata-level summary table ready for allocate_twostage_ga()
#   (2) a suggested initial c_h vector (proportional allocation)
# =============================================================================

#' Prepare a sampling frame from SSU-level census data
#'
#' Takes a dataset with one row per SSU (e.g. one row per ECD centre) and
#' aggregates it into the stratum-level summary required by
#' \code{\link{allocate_twostage_ga_joint}}. Computes stratum weights, between-
#' and within-cluster variance components, and a suggested initial cluster
#' allocation proportional to cluster count.
#'
#' @param data A data frame with one row per SSU, as returned by
#'   \code{\link{read_survey_data}}. Must contain columns named
#'   \code{strata} and \code{cluster} (standardised by \code{read_survey_data})
#'   plus the analysis variable column.
#' @param y_var Character. Name of the analysis variable column (numeric).
#'   Must match the column name retained by \code{\link{read_survey_data}}.
#' @param total_c Integer. Total number of clusters to allocate across all
#'   strata. Used to compute the suggested \code{c_h} vector.
#' @param cost_cluster Numeric vector of length equal to the number of strata,
#'   or a single value applied to all strata. Cost of sampling one cluster
#'   in each stratum. Default is \code{1} for all strata.
#' @param cost_unit Numeric vector of length equal to the number of strata,
#'   or a single value applied to all strata. Cost of sampling one SSU
#'   within a selected cluster. Default is \code{1} for all strata.
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{\code{strata}}{A data frame with one row per stratum, containing:
#'     \itemize{
#'       \item \code{stratum} — stratum identifier
#'       \item \code{C_h} — total clusters in stratum
#'       \item \code{N_h} — total SSUs in stratum
#'       \item \code{N_bar_h} — mean SSUs per cluster (\eqn{\bar{N}_h})
#'       \item \code{W_h} — stratum weight (proportion of total SSUs)
#'       \item \code{sigma2_bh} — between-cluster variance component
#'       \item \code{sigma2_wh} — within-cluster variance component
#'       \item \code{cost_cluster} — cost per sampled cluster
#'       \item \code{cost_unit} — cost per sampled SSU
#'     }
#'   }
#'   \item{\code{c_h}}{An integer vector of length \eqn{H} giving the
#'     suggested initial cluster allocation, proportional to \eqn{C_h}.}
#' }
#'
#' @details
#' The variance components follow Cochran (1977, Chapter 10):
#'
#' \deqn{\sigma^2_{bh} = \frac{1}{C_h - 1}
#'   \sum_{i=1}^{C_h} \left(\bar{y}_{hi} - \bar{\bar{y}}_h\right)^2}
#'
#' \deqn{\sigma^2_{wh} = \frac{1}{C_h}
#'   \sum_{i=1}^{C_h} \frac{1}{N_{hi} - 1}
#'   \sum_{j=1}^{N_{hi}} \left(y_{hij} - \bar{y}_{hi}\right)^2}
#'
#' where \eqn{\bar{y}_{hi}} is the cluster mean and
#' \eqn{\bar{\bar{y}}_h} is the stratum mean.
#'
#' The stratum weight is:
#' \deqn{W_h = \frac{N_h}{\sum_h N_h}}
#' where \eqn{N_h = \sum_{i=1}^{C_h} N_{hi}} is the total SSU count
#' in stratum \eqn{h}.
#'
#' The suggested initial cluster allocation is proportional to \eqn{C_h}:
#' \deqn{c_h^{(0)} = \text{round}\!\left(c_{total} \cdot
#'   \frac{C_h}{\sum_h C_h}\right)}
#' with remainders distributed to the largest strata to ensure
#' \eqn{\sum_h c_h^{(0)} = c_{total}}.
#'
#' @references
#' Cochran, W. G. (1977). \emph{Sampling Techniques} (3rd ed.). Wiley.
#'
#' @seealso \code{\link{read_survey_data}},
#'   \code{\link{allocate_twostage_ga_joint}}
#'
#' @examples
#' data("ecd_subset", package = "gaSampling")
#'
#' clean <- read_survey_data(
#'   data       = ecd_subset,
#'   stratum_id = "province",
#'   cluster_id = "id_ward",
#'   y_var      = "count_children_present"
#' )
#'
#' frame <- prepare_sampling_frame(
#'   data         = clean,
#'   y_var        = "count_children_present",
#'   total_c      = 60,
#'   cost_cluster = 1000,
#'   cost_unit    = 100
#' )
#'
#' frame$strata   # stratum summary table
#' frame$c_h      # suggested initial cluster allocation
#'
#' @export
prepare_sampling_frame <- 
  function(data,
           y_var,
           total_c,
           cost_cluster = 1,
           cost_unit    = 1) {
    
    # --- 1. Validate raw SSU-level input ---------------------------------------
    # data must have come through read_survey_data() — strata and cluster
    # are the fixed standardised column names
    .validate_ssu_input(data, y_var, total_c)
    
    # --- 2. PSU-level intermediates --------------------------------------------
    # One row per cluster: cluster mean, cluster size, stratum membership
    
    psu_level <- do.call(rbind, lapply(
      split(data, list(data[["strata"]], data[["cluster"]]), drop = TRUE),
      function(grp) {
        data.frame(
          stratum = grp[["strata"]][1],
          psu     = grp[["cluster"]][1],
          N_hi    = nrow(grp),                        # SSU count in cluster
          ybar_hi = mean(grp[[y_var]], na.rm = TRUE), # cluster mean
          stringsAsFactors = FALSE
        )
      }
    ))
    
    # --- 3. Stratum-level aggregation ------------------------------------------
    
    strata_list <- lapply(
      split(psu_level, psu_level$stratum),
      function(h) {
        
        C_h     <- nrow(h)                   # total clusters in stratum
        N_h     <- sum(h$N_hi)               # total SSUs in stratum
        N_bar_h <- mean(h$N_hi)              # mean cluster size
        
        # Grand mean of stratum (SSU-weighted)
        ybar_bar_h <- sum(h$N_hi * h$ybar_hi) / N_h
        
        # Between-cluster variance
        sigma2_bh <- if (C_h > 1) {
          sum((h$ybar_hi - ybar_bar_h)^2) / (C_h - 1)
        } else {
          NA_real_
        }
        
        # Within-cluster variance: average of per-cluster sample variances
        # We recover sum of squares from cluster mean and size using the
        # per-cluster SSU records from the original data
        # (passed through via environment — see .compute_sigma2_wh)
        sigma2_wh <- .compute_sigma2_wh(
          data    = data,
          y_var   = y_var,
          stratum = h$stratum[1],
          psu_vec = h$psu
        )
        
        data.frame(
          stratum   = h$stratum[1],
          C_h       = C_h,
          N_h       = N_h,
          N_bar_h   = N_bar_h,
          W_h       = NA_real_,       # filled after all strata computed
          sigma2_bh = sigma2_bh,
          sigma2_wh = sigma2_wh,
          stringsAsFactors = FALSE
        )
      }
    )
    
    strata <- do.call(rbind, strata_list)
    rownames(strata) <- NULL
    
    # Singleton-only strata can produce undefined variance components.
    strata$sigma2_bh[!is.finite(strata$sigma2_bh)] <- 0
    strata$sigma2_wh[!is.finite(strata$sigma2_wh)] <- 0
    
    # --- 4. Stratum weights (proportion of total SSUs) -------------------------
    strata$W_h <- strata$N_h / sum(strata$N_h)
    
    # --- 5. Cost columns -------------------------------------------------------
    H <- nrow(strata)
    
    strata$cost_cluster <- .recycle_cost(cost_cluster, H, "cost_cluster")
    strata$cost_unit    <- .recycle_cost(cost_unit,    H, "cost_unit")
    
    # --- 6. Suggested c_h: proportional to C_h ---------------------------------
    c_h <- .proportional_ch(strata$C_h, total_c)
    
    # --- 7. Print summary ------------------------------------------------------
    .print_frame_summary(strata, c_h)
    
    # --- 8. Return -------------------------------------------------------------
    list(
      strata = strata,
      c_h    = c_h
    )
  }


# =============================================================================
# Internal helpers
# =============================================================================

# --- .validate_ssu_input() ---------------------------------------------------
# Validates raw SSU-level input before any aggregation

.validate_ssu_input <- function(data, y_var, total_c) {
  
  if (!is.data.frame(data))
    stop("`data` must be a data frame.", call. = FALSE)
  
  # Expect standardised names from read_survey_data()
  missing_std <- setdiff(c("strata", "cluster"), names(data))
  if (length(missing_std) > 0)
    stop(
      "Expected standardised columns not found: ",
      paste(missing_std, collapse = ", "),
      ".\nDid you run `read_survey_data()` first?",
      call. = FALSE
    )
  
  if (!y_var %in% names(data))
    stop("Column `", y_var, "` not found in `data`.", call. = FALSE)
  
  if (!is.numeric(data[[y_var]]))
    stop("`", y_var, "` must be numeric.", call. = FALSE)
  
  if (!is.numeric(total_c) || length(total_c) != 1 || total_c < 1)
    stop("`total_c` must be a single positive integer.", call. = FALSE)
  
  total_clusters <- length(unique(data[["cluster"]]))
  if (total_c > total_clusters)
    stop("`total_c` (", total_c, ") exceeds the total number of clusters (",
         total_clusters, ") in the data.", call. = FALSE)
  
  invisible(TRUE)
}


# --- .compute_sigma2_wh() ----------------------------------------------------
# Within-cluster variance for stratum h: average of per-cluster variances
# Cochran (1977) eq. 10.7

.compute_sigma2_wh <- function(data, y_var, stratum, psu_vec) {
  
  s2_per_cluster <- sapply(psu_vec, function(psu) {
    rows <- data[["strata"]] == stratum & data[["cluster"]] == psu
    y    <- data[[y_var]][rows]
    N_hi <- sum(rows)
    if (N_hi < 2) return(NA_real_)   # variance undefined for singleton clusters
    var(y, na.rm = TRUE)              # uses 1/(N_hi - 1) denominator
  })
  
  cluster_sizes <- sapply(psu_vec, function(psu) {
    sum(data[["strata"]] == stratum & data[["cluster"]] == psu)
  })
  
  weighted.mean(s2_per_cluster, cluster_sizes, na.rm = TRUE)
}


# --- .recycle_cost() ---------------------------------------------------------
# Recycles a scalar cost value to length H, or validates a vector

.recycle_cost <- function(cost, H, name) {
  if (length(cost) == 1) {
    return(rep(cost, H))
  }
  if (length(cost) != H) {
    stop("`", name, "` must be length 1 or length ", H,
         " (one value per stratum).", call. = FALSE)
  }
  cost
}


# --- .proportional_ch() ------------------------------------------------------
# Proportional allocation of total_c clusters across strata
# Remainders distributed to largest strata to hit total_c exactly

.proportional_ch <- function(C_h, total_c) {
  
  H       <- length(C_h)
  prop    <- C_h / sum(C_h)
  c_h_raw <- prop * total_c
  c_h     <- floor(c_h_raw)
  
  # Ensure minimum of 2 per stratum (need >= 2 to estimate variance)
  c_h     <- pmax(c_h, 2L)
  
  # Distribute remainder to strata with largest fractional parts
  remainder <- total_c - sum(c_h)
  if (remainder > 0) {
    frac    <- c_h_raw - floor(c_h_raw)
    top_idx <- order(frac, decreasing = TRUE)[seq_len(remainder)]
    c_h[top_idx] <- c_h[top_idx] + 1L
  }
  
  as.integer(c_h)
}


# --- .print_frame_summary() --------------------------------------------------

.print_frame_summary <- function(strata, c_h) {
  cat("\n-- gasampling: sampling frame prepared --\n\n")
  cat(sprintf("  Strata (H)        : %d\n",   nrow(strata)))
  cat(sprintf("  Total clusters    : %d\n",   sum(strata$C_h)))
  cat(sprintf("  Total SSUs        : %d\n",   sum(strata$N_h)))
  cat(sprintf("  Suggested c_h     : %s  (total = %d)\n",
              paste(c_h, collapse = " "), sum(c_h)))
  cat("\n  Stratum summary:\n")
  cat(sprintf("  %-20s %6s %8s %8s %10s %10s\n",
              "Stratum", "C_h", "N_h", "W_h", "sigma2_bh", "sigma2_wh"))
  cat("  ", strrep("-", 78), "\n", sep = "")
  for (i in seq_len(nrow(strata))) {
    cat(sprintf("  %-20s %6d %8d %8.4f %10.4f %10.4f\n",
                strata$stratum[i],
                strata$C_h[i],
                strata$N_h[i],
                strata$W_h[i],
                strata$sigma2_bh[i],
                strata$sigma2_wh[i]))
  }
  cat("\n")
}