# =============================================================================
# gasampling: read_survey_data()
# Place in: R/read_survey_data.R
#
# Takes a user's data frame with arbitrary column names and returns a
# standardised data frame with gasampling convention column names,
# ready for prepare_sampling_frame().
#
# Standard output column names:
#   strata  — stratum identifier
#   cluster — PSU/cluster identifier
#   y_var   — analysis variable (retains original name)
# =============================================================================

#' Read and standardise survey data for use with gasampling
#'
#' Takes a data frame with arbitrary column names and returns a cleaned,
#' standardised data frame using \code{gasampling} column name conventions.
#' This is the entry point for all user data — the output feeds directly
#' into \code{\link{prepare_sampling_frame}}.
#'
#' @param data A data frame with one row per SSU. Can have any column names.
#' @param stratum_id Character. Name of the column in \code{data} that
#'   identifies strata (e.g. \code{"province"}, \code{"region"}).
#' @param cluster_id Character. Name of the column in \code{data} that
#'   identifies clusters/PSUs (e.g. \code{"district"}, \code{"school"}).
#' @param y_var Character. Name of the analysis variable column (numeric).
#'   This column retains its original name in the output.
#' @param na_action Character. How to handle missing values in \code{y_var}.
#'   One of:
#'   \describe{
#'     \item{\code{"warn"}}{Drop rows with missing \code{y_var} and issue a
#'       warning with the count. Default.}
#'     \item{\code{"stop"}}{Stop with an error if any missing values found.}
#'     \item{\code{"keep"}}{Keep missing values as-is. Downstream functions
#'       may handle or flag them.}
#'   }
#'
#' @return A data frame with the same rows as \code{data} (minus any dropped
#'   missing rows) and standardised column names:
#'   \describe{
#'     \item{\code{strata}}{Stratum identifier, coerced to character,
#'       whitespace stripped.}
#'     \item{\code{cluster}}{Cluster identifier, coerced to character,
#'       whitespace stripped.}
#'     \item{\code{<y_var>}}{Analysis variable, coerced to numeric.
#'       Retains the original column name supplied by the user.}
#'   }
#'   Any additional columns in \code{data} beyond the three mapped ones
#'   are dropped.
#'
#' @details
#' Cleaning steps applied in order:
#' \enumerate{
#'   \item Checks that all named columns exist in \code{data}.
#'   \item Coerces \code{stratum_id} and \code{cluster_id} to character
#'     and strips leading/trailing whitespace.
#'   \item Checks that clusters nest cleanly within strata (no cluster
#'     appearing in more than one stratum).
#'   \item Coerces \code{y_var} to numeric, warning if coercion introduces
#'     new \code{NA}s.
#'   \item Handles missing values in \code{y_var} according to
#'     \code{na_action}.
#'   \item Checks for negative values in \code{y_var} and warns if found.
#'   \item Checks for singleton clusters (clusters with only one SSU) and
#'     warns if found, as within-cluster variance is undefined for these.
#'   \item Renames columns to \code{gasampling} standard names.
#' }
#'
#' @seealso \code{\link{prepare_sampling_frame}}
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
#' head(clean)
#'
#' @export
read_survey_data <- function(data,
                             stratum_id,
                             cluster_id,
                             y_var,
                             na_action = c("warn", "stop", "keep")) {
  
  na_action <- match.arg(na_action)
  
  # --- 1. Input is a data frame ----------------------------------------------
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  
  # --- 2. Named columns exist ------------------------------------------------
  required <- c(stratum_id, cluster_id, y_var)
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0) {
    stop(
      "The following columns were not found in `data`:\n  ",
      paste(missing_cols, collapse = ", "),
      "\nAvailable columns:\n  ",
      paste(names(data), collapse = ", "),
      call. = FALSE
    )
  }
  
  # --- 3. Subset to only the columns we need ---------------------------------
  out <- data[, c(stratum_id, cluster_id, y_var), drop = FALSE]
  
  # --- 4. Clean stratum and cluster IDs --------------------------------------
  out[[stratum_id]] <- trimws(as.character(out[[stratum_id]]))
  out[[cluster_id]] <- trimws(as.character(out[[cluster_id]]))
  
  # Check for empty strings after trimming
  if (any(out[[stratum_id]] == "")) {
    stop("`", stratum_id, "` contains empty strings after whitespace trimming.",
         call. = FALSE)
  }
  if (any(out[[cluster_id]] == "")) {
    stop("`", cluster_id, "` contains empty strings after whitespace trimming.",
         call. = FALSE)
  }
  
  # --- 5. Clusters nest within strata ----------------------------------------
  cluster_strata <- tapply(
    out[[stratum_id]],
    out[[cluster_id]],
    function(x) length(unique(x))
  )
  crossed <- names(cluster_strata[cluster_strata > 1])
  if (length(crossed) > 0) {
    stop(
      "The following clusters appear in more than one stratum:\n  ",
      paste(head(crossed, 5), collapse = ", "),
      if (length(crossed) > 5)
        paste0(" ... and ", length(crossed) - 5, " more"),
      "\nEach cluster must belong to exactly one stratum.",
      call. = FALSE
    )
  }
  # --- 6. Coerce y_var to numeric --------------------------------------------
  y_original <- out[[y_var]]
  y_numeric  <- suppressWarnings(as.numeric(as.character(out[[y_var]])))
  
  # Safely check if coercion introduced new NAs
  coercion_nas <- is.na(y_numeric) & !is.na(y_original)
  coercion_nas <- coercion_nas & !is.na(coercion_nas) # Remove potential NAs from logical check
  
  out[[y_var]] <- y_numeric
  
  if (isTRUE(any(coercion_nas))) {
    warning(
      sum(coercion_nas), " value(s) in `", y_var,
      "` could not be coerced to numeric and were set to NA.\n",
      "Check for non-numeric entries in this column.",
      call. = FALSE
    )
  }
  # --- 7. Handle missing values in y_var -------------------------------------
  n_missing <- sum(is.na(out[[y_var]]))
  
  if (n_missing > 0) {
    msg <- paste0(
      n_missing, " row(s) (", round(100 * n_missing / nrow(out), 1),
      "%) have missing values in `", y_var, "`."
    )
    if (na_action == "stop") {
      stop(msg, " Set `na_action = 'warn'` to drop them or ",
           "`na_action = 'keep'` to retain them.", call. = FALSE)
    } else if (na_action == "warn") {
      warning(msg, " These rows have been dropped.", call. = FALSE)
      out <- out[!is.na(out[[y_var]]), ]
    }
    # na_action == "keep": do nothing
  }
  
  # --- 8. Check for negative values ------------------------------------------
  n_negative <- sum(out[[y_var]] < 0, na.rm = TRUE)
  if (n_negative > 0) {
    warning(
      n_negative, " row(s) have negative values in `", y_var, "`. ",
      "Check whether these are valid for your analysis variable.",
      call. = FALSE
    )
  }
  
  # --- 9. Singleton cluster warning ------------------------------------------
  cluster_sizes <- table(out[[cluster_id]])
  singletons <- names(cluster_sizes[cluster_sizes == 1])
  if (length(singletons) > 0) {
    warning(
      length(singletons), " cluster(s) contain only one SSU:\n  ",
      paste(head(singletons, 5), collapse = ", "),
      if (length(singletons) > 5)
        paste0(" ... and ", length(singletons) - 5, " more"),
      "\nWithin-cluster variance is undefined for singleton clusters. ",
      "Consider removing them or merging with adjacent clusters.",
      call. = FALSE
    )
  }
  
  # --- 10. Rename to gasampling standard names --------------------------------
  names(out)[names(out) == stratum_id] <- "strata"
  names(out)[names(out) == cluster_id] <- "cluster"
  # y_var retains its original name
  
  # --- 11. Reset row names and print summary ---------------------------------
  rownames(out) <- NULL
  .print_read_summary(out, y_var)
  
  out
}


# =============================================================================
# Internal helper
# =============================================================================

.print_read_summary <- function(data, y_var) {
  cat("\n-- gasampling: survey data loaded --\n\n")
  cat(sprintf("  Rows (SSUs)    : %d\n",  nrow(data)))
  cat(sprintf("  Strata         : %d\n",  length(unique(data$strata))))
  cat(sprintf("  Clusters       : %d\n",  length(unique(data$cluster))))
  cat(sprintf("  Analysis var   : %s\n",  y_var))
  cat(sprintf("  y range        : %.2f -- %.2f\n",
              min(data[[y_var]], na.rm = TRUE),
              max(data[[y_var]], na.rm = TRUE)))
  cat(sprintf("  Missing in y   : %d\n",  sum(is.na(data[[y_var]]))))
  cat("\n  Standardised column names: strata | cluster |", y_var, "\n\n")
}