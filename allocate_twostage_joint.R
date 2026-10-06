#' Joint two-stage genetic algorithm allocation function
#'
#' Optimise cluster and unit allocations across strata
#'
#' Use a genetic algorithm to jointly search for cluster counts and
#' unit counts per cluster that minimise the variance of the overall
#' mean, subject to allocation bounds, a minimum sample size and
#' an optional budget range. Repair candidate allocations to obtain
#' integer counts and attempt to satisfy the constraints.
#'
#' @param strata A data frame with one row per stratum containing
#'   stratum (stratum identifier), C_h (available cluster count),
#'   W_h (population weight), sigma2_bh (between-cluster variance),
#'   sigma2_wh (within-cluster variance), cost_cluster (cost per
#'   sampled cluster) and cost_unit (cost per sampled secondary unit).
#'   Each stratum must have at least two available clusters.
#' @param total_n A single positive integer giving the minimum required
#'   total number of sampled secondary units.
#' @param max_cost A single positive numeric value giving the reference
#'   budget. Supply NULL to disable budget constraints.
#' @param popSize Positive integer giving the genetic algorithm
#'   population size. Defaults to 100.
#' @param maxiter Positive integer giving the maximum number of genetic
#'   algorithm iterations. Defaults to 1000.
#' @param run Positive integer giving the number of consecutive
#'   iterations without improvement in the best fitness before
#'   stopping. Defaults to 50.
#' @param seed Optional integer random seed passed to set.seed before
#'   the search. Defaults to NULL, leaving the current random-number
#'   state unchanged before the search begins.
#' @param suggestion Optional numeric vector of length 2 * nrow(strata)
#'   supplying one starting allocation. The first half contains cluster
#'   counts and the second half contains unit counts per cluster,
#'   both in the row order of strata. If NULL, a starting allocation
#'   is constructed and repaired automatically.
#' @param monitor Monitoring argument passed to GA::ga. Defaults to
#'   NULL. A monitoring function can be supplied to report search
#'   progress.
#'
#' @return A named list containing:
#'   \item{feasible}{Logical value indicating whether the returned
#'     allocation satisfies the checked constraints.}
#'   \item{fallback_used}{Logical value indicating whether a feasible
#'     repaired version of the supplied suggestion replaced an
#'     infeasible repaired genetic algorithm solution.}
#'   \item{variance}{Total variance of the overall mean for the
#'     returned allocation.}
#'   \item{between_variance}{Total weighted between-cluster contribution
#'     to the variance.}
#'   \item{within_variance}{Total weighted within-cluster contribution
#'     to the variance.}
#'   \item{total_n}{Achieved total sample size, which may exceed the
#'     requested minimum.}
#'   \item{total_cost}{Total sampling cost of the returned allocation.}
#'   \item{allocation}{Data frame containing stratum, sampled cluster
#'     counts c_h, sampled unit counts per cluster k_h, stratum sample
#'     sizes n_h, and weighted variance contributions var_between,
#'     var_within, and var_total.}
#'   \item{constraint_check}{Diagnostic list returned by
#'     check_constraints for the final allocation.}
#'   \item{repair}{Repair results and diagnostics for the returned
#'     allocation.}
#'   \item{ga_object}{The fitted GA object returned by GA::ga.}
#'
#' @details
#'   The search uses real-valued candidates, which are rounded and
#'   repaired during fitness evaluation. Cluster counts are bounded
#'   between 2 and strata$C_h. Unit counts per cluster are bounded
#'   between 1 and total_n.
#'
#'   The upper unit bound does not represent observed cluster capacity.
#'   The function does not check whether individual clusters contain
#'   enough secondary units for the returned allocation.
#'
#'   Fitness is the negative variance for feasible repaired allocations,
#'   with a fixed penalty for infeasible allocations. Variance is
#'   calculated without finite population corrections.
#'
#'   With the current helper defaults, budget constraints require
#'   spending between 95 and 102 percent of max_cost. These proportions
#'   are not exposed as arguments to this function.
#'
#'   The first solution stored in the fitted GA object is selected
#'   and repaired. If it remains infeasible and a suggestion was
#'   supplied, the function attempts to use a feasible repaired version
#'   of that suggestion instead.
#'
#'   The search is stochastic and does not guarantee a globally optimal
#'   or feasible allocation. An infeasible result can still be returned.
#'   Check feasible and constraint_check before using it. The returned
#'   allocation may differ from the solution stored in ga_object
#'   because of final repair or fallback.
#'
#'   The function requires the GA package and the helper functions
#'   fitness_ckh_constrained, repair_joint, check_constraints,
#'   cost_twostage and variance_twostage.
#'   Supplying seed resets the session's random-number state and
#'   running the search advances it. Input validation is limited.
#'
#' @examples
#' strata <- data.frame(
#'   stratum = c("A", "B"),
#'   C_h = c(20, 25),
#'   W_h = c(0.4, 0.6),
#'   sigma2_bh = c(2, 3),
#'   sigma2_wh = c(6, 9),
#'   cost_cluster = c(100, 120),
#'   cost_unit = c(5, 6)
#' )
#'
#' result <- allocate_twostage_ga_joint(
#'   strata = strata,
#'   total_n = 150,
#'   max_cost = 3900,
#'   popSize = 50,
#'   maxiter = 100,
#'   run = 20,
#'   seed = 123,
#'   suggestion = c(10, 15, 5, 8)
#' )
#'
#' result$feasible
#' result$allocation
#' result$constraint_check
allocate_twostage_ga_joint <- function(strata,
                                       total_n,
                                       max_cost,
                                       popSize = 100,
                                       maxiter = 1000,
                                       run = 50,
                                       seed = NULL,
                                       suggestion = NULL,
                                       monitor = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  L <- nrow(strata)
  
  if (!"C_h" %in% names(strata)) {
    stop("strata must contain C_h, the number of available clusters per stratum.")
  }
  
  lower_c <- rep(2L, L)
  lower_k <- rep(1L, L)
  upper_c <- strata$C_h
  upper_k <- rep(total_n, L)
  
  max_total_n <- sum(upper_c * upper_k)
  
  if (total_n > max_total_n) {
    stop(
      "The requested minimum sample size exceeds the maximum ",
      "possible sample size under the allocation bounds: ",
      max_total_n, "."
    )
  }
  
  lower <- c(
    lower_c,
    lower_k
  )
  
  upper <- c(
    upper_c,
    upper_k
  )
  
  if (is.null(suggestion)) {
    # if no suggestion, then propose proportional for each candidate
    starting_c_h <- pmax(lower_c, pmin(upper_c, round(total_n / L / 2)))
    starting_k_h <- pmax(lower_k, round(total_n / sum(starting_c_h)))
    
    starting <- repair_joint(
      c_h = starting_c_h,
      k_h = starting_k_h,
      total_n = total_n,
      lower_c = lower_c,
      upper_c = upper_c,
      lower_k = lower_k,
      upper_k = upper_k,
      max_cost = max_cost,
      cost_cluster = strata$cost_cluster,
      cost_unit = strata$cost_unit
    )
    
    suggestions <- matrix(
      c(starting$c_h, starting$k_h),
      nrow = 1
    )
    
  } else {
    
    suggestions <- matrix(
      suggestion,
      nrow = 1
    )
  }
  
  ga_fit <- GA::ga(
    type = "real-valued",
    fitness = function(x) {
      fitness_ckh_constrained(
        x = x,
        strata = strata,
        total_n = total_n,
        max_cost = max_cost,
        min_c_h = 2L,
        min_k_h = 1L,
        max_c_h = upper_c,
        max_k_h = upper_k
      )
    },
    lower = lower,
    upper = upper,
    suggestions = suggestions,
    popSize = popSize,
    maxiter = maxiter,
    run = run,
    monitor = monitor
  )
  
  best_x <- ga_fit@solution[1, ]
  
  best_c_h <- round(best_x[1:L])
  best_k_h <- round(best_x[(L + 1):(2 * L)])
  
  # repair final solution
  repaired <- repair_joint(
    c_h = best_c_h,
    k_h = best_k_h,
    total_n = total_n,
    lower_c = lower_c,
    upper_c = upper_c,
    lower_k = lower_k,
    upper_k = upper_k,
    max_cost = max_cost,
    cost_cluster = strata$cost_cluster,
    cost_unit = strata$cost_unit
  )
  
  repaired_check <- check_constraints(
    c_h = repaired$c_h,
    k_h = repaired$k_h,
    total_n = total_n,
    max_cost = max_cost,
    cost_cluster = strata$cost_cluster,
    cost_unit = strata$cost_unit,
    min_clusters = 2L,
    min_units = 1L,
    max_clusters = upper_c,
    max_units = upper_k
  )
  
  fallback_used <- FALSE # Track if a fallback was used
  
  # Fallback condition
  if (!repaired_check$feasible && !is.null(suggestion)) {
    suggestion_repaired <- repair_joint(
      c_h = suggestion[1:L],
      k_h = suggestion[(L + 1):(2 * L)],
      total_n = total_n,
      lower_c = lower_c,
      upper_c = upper_c,
      lower_k = lower_k,
      upper_k = upper_k,
      max_cost = max_cost,
      cost_cluster = strata$cost_cluster,
      cost_unit = strata$cost_unit
    )
    
    suggestion_check <- check_constraints(
      c_h = suggestion_repaired$c_h,
      k_h = suggestion_repaired$k_h,
      total_n = total_n,
      max_cost = max_cost,
      cost_cluster = strata$cost_cluster,
      cost_unit = strata$cost_unit,
      min_clusters = 2L,
      min_units = 1L,
      max_clusters = upper_c,
      max_units = upper_k
    )
    
    if (suggestion_check$feasible) {
      repaired <- suggestion_repaired
      fallback_used <- TRUE # When the starting suggestion replaces the solution
    }
  }
  
  best_c_h <- repaired$c_h
  best_k_h <- repaired$k_h

  allocation <- data.frame(
    stratum = strata$stratum,
    c_h = best_c_h,
    k_h = best_k_h,
    n_h = best_c_h * best_k_h
  ) 
  

  var_between_h <- (strata$W_h^2 * strata$sigma2_bh) / best_c_h
  var_within_h  <- (strata$W_h^2 * strata$sigma2_wh) / (best_c_h * best_k_h)
  var_total_h   <- var_between_h + var_within_h
  
  allocation$var_between <- var_between_h
  allocation$var_within  <- var_within_h
  allocation$var_total   <- var_total_h
  
  # Global sums for the summary statistics
  total_variance <- sum(var_total_h)
  total_between  <- sum(var_between_h)
  total_within   <- sum(var_within_h)
  
  total_cost <- cost_twostage(
    c_h = best_c_h,
    k_h = best_k_h,
    cost_cluster = strata$cost_cluster,
    cost_unit = strata$cost_unit
  ) 
  
  constraint_check <- check_constraints(
    c_h = best_c_h,
    k_h = best_k_h,
    total_n = total_n,
    max_cost = max_cost,
    cost_cluster = strata$cost_cluster,
    cost_unit = strata$cost_unit,
    min_clusters = 2L,
    min_units = 1L,
    max_clusters = upper_c,
    max_units = upper_k
  )
  
  is_feasible <- constraint_check$feasible && repaired$repaired_exactly
  
  list(
    feasible         = is_feasible,
    fallback_used    = fallback_used,
    variance         = total_variance,
    between_variance = total_between,    
    within_variance  = total_within,     
    total_n          = sum(best_c_h * best_k_h),  
    total_cost       = total_cost,  
    allocation       = allocation,
    constraint_check = constraint_check,  
    repair           = repaired,  
    ga_object        = ga_fit  
  )
}