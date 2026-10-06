#' Two-stage allocation constraint checking function
#'
#' Check whether a two-stage allocation satisfies sampling constraints
#'
#' Check integer cluster and unit counts, allocation bounds, an optional
#' minimum total sample size and an optional budget range.
#' Return overall feasibility and separate indicators for each check.
#'
#' @param c_h Numeric vector giving the number of primary sampling units
#'   (clusters) sampled in each stratum.
#' @param k_h Numeric vector giving the number of secondary sampling
#'   units sampled per cluster in each stratum.
#' @param total_n A single non-negative integer giving the minimum
#'   required total sample size. Defaults to NULL, which disables
#'   the sample size check.
#' @param max_cost A single positive numeric value giving the reference
#'   budget. Defaults to NULL, which disables the budget check.
#' @param cost_cluster Numeric vector giving the cost per sampled
#'   cluster in each stratum, or a single value applied to all strata.
#'   Values should be finite and non-negative. Required when max_cost
#'   is supplied.
#' @param cost_unit Numeric vector giving the cost per sampled secondary
#'   unit in each stratum, or a single value applied to all strata.
#'   Values should be finite and non-negative. Required when max_cost
#'   is supplied.
#' @param min_clusters A single integer or a vector of integers giving
#'   the minimum cluster count per stratum. Defaults to 2.
#' @param min_units A single integer or a vector of integers giving
#'   the minimum number of sampled units per cluster in each stratum.
#'   Defaults to 1.
#' @param max_clusters A single integer or a vector of integers giving
#'   the maximum cluster count per stratum. Defaults to NULL, which
#'   disables the upper cluster bound check.
#' @param max_units A single integer or a vector of integers giving
#'   the maximum number of sampled units per cluster in each stratum.
#'   Defaults to NULL, which disables the upper unit bound check.
#' @param min_budget_use Numeric proportion giving the minimum required
#'   budget use. Defaults to 0.95, requiring spending of at least
#'   95 percent of max_cost.
#' @param budget_tol Non-negative numeric proportion allowing spending
#'   above max_cost. Defaults to 0.02, allowing spending up to
#'   102 percent of max_cost.
#'
#' @return A named list containing:
#'   \item{feasible}{Logical value indicating whether all enabled
#'     constraints are satisfied.}
#'   \item{integer_ok}{Logical value indicating whether all cluster
#'     and unit counts are integer-valued.}
#'   \item{cluster_ok}{Logical value indicating whether cluster counts
#'     satisfy the lower and any supplied upper bounds.}
#'   \item{unit_ok}{Logical value indicating whether unit counts per
#'     cluster satisfy the lower and any supplied upper bounds.}
#'   \item{sample_size_ok}{Logical value indicating whether the minimum
#'     total sample size is met or exceeded. TRUE when total_n is NULL.}
#'   \item{cost_ok}{Logical value indicating whether spending falls
#'     within the allowed budget range. TRUE when max_cost is NULL.}
#'   \item{total_n}{The achieved total sample size, sum(c_h * k_h).}
#'   \item{total_cost}{The total sampling cost, or NA when max_cost
#'     is NULL.}
#'   \item{n_h}{Numeric vector of stratum sample sizes, c_h * k_h.}
#'
#' @details
#'   Allocation vectors should be finite and have the same length,
#'   with entries ordered consistently by stratum. Bounds and costs
#'   should each have one entry per stratum or a single common value.
#'
#'   The sample size requirement is a minimum, rather than an exact
#'   target. When budget checking is enabled, spending must lie between
#'   max_cost * min_budget_use and max_cost * (1 + budget_tol),
#'   including both endpoints.
#'
#'   Integer checks use exact equality with rounded values, without
#'   a numerical tolerance. The function does not round, repair, or
#'   otherwise modify the supplied allocation.
#'
#'   The function requires cost_twostage when max_cost is supplied.
#'   It checks that required cost inputs are present, but does not
#'   validate input lengths, missing values, or bound consistency.
#'   Invalid inputs may produce missing indicators or errors.
#'
#' @examples
#' check_constraints(
#'   c_h = c(10, 15),
#'   k_h = c(5, 8),
#'   total_n = 150,
#'   max_cost = 3900,
#'   cost_cluster = c(100, 120),
#'   cost_unit = c(5, 6),
#'   min_clusters = 2,
#'   min_units = 1,
#'   max_clusters = c(20, 25),
#'   max_units = c(10, 10)
#' )
check_constraints <- function(c_h,
                              k_h,
                              total_n = NULL,
                              max_cost = NULL,
                              cost_cluster = NULL,
                              cost_unit = NULL,
                              min_clusters = 2,
                              min_units = 1,
                              max_clusters = NULL,
                              max_units = NULL,
                              min_budget_use = 0.95,
                              budget_tol = 0.02) {
  
  # Compute stratum sample sizes
  n_h <- c_h * k_h
  
  # Check integer values
  integer_ok <- all(c_h == round(c_h)) && all(k_h == round(k_h))
  
  # Check lower bounds
  cluster_ok <- all(c_h >= min_clusters)
  unit_ok <- all(k_h >= min_units)
  
  # Check upper bounds
  if (!is.null(max_clusters)) {
    cluster_ok <- cluster_ok && all(c_h <= max_clusters)
  }
  
  if (!is.null(max_units)) {
    unit_ok <- unit_ok && all(k_h <= max_units)
  }
  
  # Check total sample size constraint
  if (!is.null(total_n)) {
    sample_size_gap <- sum(n_h) - total_n
    sample_size_ok <- sample_size_gap >= 0
  } else {
    sample_size_gap <- NA_integer_
    sample_size_ok <- TRUE
  }
  
  # Check cost constraint (FORCE A STRICTER LB WITH A SMALL UB tolerance)
  if (!is.null(max_cost)) {
    
    if (is.null(cost_cluster) || is.null(cost_unit)) {
      stop("cost_cluster and cost_unit must be supplied when max_cost is used.")
    }
    
    total_cost <- cost_twostage(c_h, k_h, cost_cluster, cost_unit)
    
    lower_cost_ok <- total_cost >= max_cost * min_budget_use
    upper_cost_ok <- total_cost <= max_cost * (1 + budget_tol)
    
    cost_ok <- lower_cost_ok && upper_cost_ok
    
  } else { # default if left null
    total_cost <- NA
    cost_ok <- TRUE
  }
  
  feasible <- integer_ok && # must all be satisfied to be a feasible solution
    cluster_ok &&
    unit_ok &&
    sample_size_ok &&
    cost_ok
  
  list(
    feasible = feasible,
    integer_ok = integer_ok,
    cluster_ok = cluster_ok,
    unit_ok = unit_ok,
    sample_size_ok = sample_size_ok,
    cost_ok = cost_ok,
    total_n = sum(n_h),
    total_cost = total_cost,
    n_h = n_h
  )
}