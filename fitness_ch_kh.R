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
