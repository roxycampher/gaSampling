#' Constrained joint allocation fitness function
#'
#' Evaluate the fitness of a two-stage sampling allocation
#'
#' Round and bound candidate cluster and unit counts, attempt to repair
#' the allocation and check sampling and budget constraints.
#' Return the negative variance for a feasible repaired allocation
#' or a fixed penalty for an infeasible allocation.
#'
#' @param x Numeric vector of length 2 * nrow(strata). The first half
#'   contains candidate cluster counts by stratum and the second half
#'   contains candidate unit counts per cluster in the same stratum order.
#' @param strata A data frame with one row per stratum containing
#'   W_h (population weight), sigma2_bh (between-cluster variance),
#'   sigma2_wh (within-cluster variance), cost_cluster (cost per sampled
#'   cluster) and cost_unit (cost per sampled secondary unit).
#'   The column C_h (available cluster count) is also required when
#'   max_c_h is NULL.
#' @param total_n A single positive integer giving the minimum required
#'   total sample size.
#' @param max_cost A single positive numeric value giving the reference
#'   budget. Supply NULL to disable budget constraints.
#' @param min_c_h A single positive integer giving the minimum number
#'   of sampled clusters in each stratum. Defaults to 2.
#' @param min_k_h A single positive integer giving the minimum number
#'   of sampled units per cluster in each stratum. Defaults to 1.
#' @param max_c_h Numeric vector giving the maximum number of sampled
#'   clusters in each stratum. Defaults to strata$C_h when NULL.
#' @param max_k_h Numeric vector giving the maximum number of sampled
#'   units per cluster in each stratum. Defaults to total_n in every
#'   stratum when NULL.
#' @param min_budget_use Numeric proportion giving the minimum required
#'   budget use. Defaults to 0.95, requiring spending of at least
#'   95 percent of max_cost.
#' @param budget_tol Non-negative numeric proportion allowing spending
#'   above max_cost. Defaults to 0.02, allowing spending up to
#'   102 percent of max_cost.
#'
#' @return A single numeric fitness score. Returns the negative variance
#'   of the repaired allocation when all checked constraints are
#'   satisfied, or -1e10 when the repaired allocation is infeasible.
#'
#' @details
#'   Candidate values are rounded to integers and restricted to their
#'   allocation bounds before repair. Upper-bound vectors must have
#'   one entry per stratum, ordered consistently with strata and x.
#'   Bounds should be integer-valued and internally consistent.
#'
#'   Feasibility requires integer allocations within the specified
#'   bounds, a total sample size at least as large as total_n and,
#'   when enabled, spending between max_cost * min_budget_use and
#'   max_cost * (1 + budget_tol).
#'
#'   The fitness is calculated from the repaired allocation, which may
#'   differ from the supplied candidate. A larger fitness score indicates
#'   a smaller variance among feasible allocations. The function returns
#'   only the score, not the repaired allocation or constraint diagnostics.
#'
#'   Repair is heuristic and may fail to find a feasible allocation.
#'   The function requires repair_joint, check_constraints, cost_twostage,
#'   and variance_twostage to be available. It does not fully validate
#'   inputs, and invalid inputs may produce an error rather than a penalty.
#'
#' @examples
#' strata <- data.frame(
#'   C_h = c(20, 25),
#'   W_h = c(0.4, 0.6),
#'   sigma2_bh = c(2, 3),
#'   sigma2_wh = c(6, 9),
#'   cost_cluster = c(100, 120),
#'   cost_unit = c(5, 6)
#' )
#'
#' fitness_ckh_constrained(
#'   x = c(10, 15, 5, 8),
#'   strata = strata,
#'   total_n = 150,
#'   max_cost = 3900,
#'   max_k_h = c(10, 10)
#' )

fitness_ckh_constrained <- function(x,
                                    strata,
                                    total_n,
                                    max_cost,
                                    min_c_h = 2L,
                                    min_k_h = 1L,
                                    max_c_h = NULL,
                                    max_k_h = NULL,
                                    min_budget_use = 0.95, # default settings
                                    budget_tol = 0.02) {
  
  L <- nrow(strata)
  
  c_h <- round(x[1:L])
  k_h <- round(x[(L + 1):(2 * L)])
  
  if (is.null(max_c_h)) max_c_h <- strata$C_h
  if (is.null(max_k_h)) max_k_h <- rep(total_n, L)
  
  c_h <- pmax(c_h, min_c_h)
  k_h <- pmax(k_h, min_k_h)
  c_h <- pmin(c_h, max_c_h)
  k_h <- pmin(k_h, max_k_h)
  
  repaired <- repair_joint(
    c_h = c_h,
    k_h = k_h,
    total_n = total_n, 
    lower_c = rep(min_c_h, L),
    upper_c = max_c_h,
    lower_k = rep(min_k_h, L),
    upper_k = max_k_h,
    max_cost = max_cost,
    cost_cluster = strata$cost_cluster,
    cost_unit = strata$cost_unit,
    min_budget_use = min_budget_use,
    budget_tol = budget_tol
  )
  
  c_h <- repaired$c_h
  k_h <- repaired$k_h
  
  constraint_check <- check_constraints(
    k_h = k_h,
    c_h = c_h,
    total_n = total_n,
    max_cost = max_cost,
    cost_cluster = strata$cost_cluster,
    cost_unit = strata$cost_unit,
    min_clusters = min_c_h,
    min_units = min_k_h,
    max_clusters = max_c_h,
    max_units = max_k_h,
    min_budget_use = min_budget_use,
    budget_tol = budget_tol
  )

  if (!constraint_check$feasible) {
    return(-1e10)
  }
  
  variance <- variance_twostage(
    c_h = c_h,
    k_h = k_h,
    W_h = strata$W_h,
    sigma2_bh = strata$sigma2_bh,
    sigma2_wh = strata$sigma2_wh
  )
  
  -variance
}
