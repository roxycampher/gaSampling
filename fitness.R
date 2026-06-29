fitness_kh_constrained <- function(x, strata, total_n, max_cost, c_h) {
  
  # x = c(k_1, k_2, ..., k_L)
  k_h <- repair_kh(
    x = x,
    c_h = c_h,
    total_n = total_n,
    lower = rep(1L, length(c_h)),
    upper = floor(total_n / c_h)
  )
  
  constraints <- check_constraints(
    c_h = c_h,
    k_h = k_h,
    total_n = total_n,
    max_cost = max_cost,
    cost_cluster = strata$cost_cluster,
    cost_unit = strata$cost_unit
  )
  
  if (!constraints$feasible) {
    return(-1e15)
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