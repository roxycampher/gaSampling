fitness_ckh_constrained <- function(x,
                                    strata,
                                    total_n,
                                    max_cost,
                                    min_c_h = 1L,
                                    min_k_h = 1L,
                                    max_c_h = NULL,
                                    max_k_h = NULL,
                                    sample_shortfall_penalty = 1e5,
                                    cost_over_penalty = 1e5,
                                    budget_underuse_penalty = 1e3,
                                    min_budget_use = 0.95,
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
  
  # KEEP total_n passed here since your custom repair logic now handles >= total_n!
  repaired <- repair_ckh(
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
  
  # KEEP total_n here too so check_constraints evaluates against your minimum boundary
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
  

#  if (!constraint_check$integer_ok ||
#      !constraint_check$cluster_ok ||
#      !constraint_check$unit_ok) {

  # If repair function changed how it flags 'repaired_exactly' (e.g. meaning 'successfully >= total_n'),
  # keep it! Otherwise, if check_constraints returns feasible = TRUE, we trust it.
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
  
  # sample_shortfall <- max(0, total_n - constraint_check$total_n)
  # upper_cost_bound <- max_cost * (1 + budget_tol)
  # lower_cost_bound <- max_cost * min_budget_use
  # over_cost <- max(0, constraint_check$total_cost - upper_cost_bound)
  # budget_underuse <- max(0, lower_cost_bound - constraint_check$total_cost)
  
  -variance
}
