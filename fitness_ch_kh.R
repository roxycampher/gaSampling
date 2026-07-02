fitness_ckh_constrained <- function(x,
                                    strata,
                                    total_n,
                                    max_cost,
                                    min_c_h = 1L,
                                    min_k_h = 1L) {
  
  L <- nrow(strata)
  
  c_h <- round(x[1:L])
  k_h <- round(x[(L + 1):(2 * L)])
  
  c_h <- pmax(c_h, min_c_h)
  k_h <- pmax(k_h, min_k_h)
  
  # repair c_h and k_h jointly
  repaired <- repair_ckh(
    c_h = c_h,
    k_h = k_h,
    total_n = total_n,
    lower_c = rep(min_c_h, L),
    lower_k = rep(min_k_h, L)
  )
  
  c_h <- repaired$c_h
  k_h <- repaired$k_h
  
  constraint_check <- check_constraints(
    k_h = k_h,
    c_h = c_h,
    total_n = total_n,
    max_cost = max_cost,
    cost_cluster = strata$cost_cluster,
    cost_unit = strata$cost_unit
  )
  
  if (!constraint_check$valid) {
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