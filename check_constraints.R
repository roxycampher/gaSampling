# so that users can debug their solution, where it goes wrong etc
check_constraints <- function(c_h,
                              k_h,
                              total_n = NULL,
                              max_cost = NULL,
                              cost_cluster = NULL,
                              cost_unit = NULL,
                              min_clusters = 1,
                              min_units = 1,
                              max_clusters = NULL,
                              max_units = NULL,
                              min_budget_use = NULL,
                              budget_tol = NULL) {
  
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
  
  # Check total sample size constraint (LOWER BOUND!!)
  if (!is.null(total_n)) {
    sample_size_ok <- sum(n_h) >= total_n
  } else {
    sample_size_ok <- TRUE # if no total sample size was set
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