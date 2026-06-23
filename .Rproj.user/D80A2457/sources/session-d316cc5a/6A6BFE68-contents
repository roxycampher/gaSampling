fitness_constrained <- function(x, strata, total_n, max_cost,
                                min_clusters = 1,
                                min_units = 1) {
  
  L <- nrow(strata)
  
  # Force integer candidate solution (migth have to change)
  x <- round(x)
  
  c_h <- x[1:L]
  k_h <- x[(L + 1):(2 * L)]
  n_h <- c_h * k_h
  
  # Constraint 1: positive lower bounds
  if (any(c_h < min_clusters) || any(k_h < min_units)) {
    return(-1e15)
  }
  
  # Constraint 2: total sample size exactly equals total_n
  if (sum(n_h) != total_n) {
    return(-1e15)
  }
  
  # Constraint 3: cost must not exceed budget
  total_cost <- cost_twostage(
    c_h = c_h,
    k_h = k_h,
    cost_cluster = strata$cost_cluster,
    cost_unit = strata$cost_unit
  )
  
  if (total_cost > max_cost) {
    return(-1e15)
  }
  
  # Objective: minimise variance, so maximise negative variance
  variance <- variance_twostage(
    c_h = c_h,
    k_h = k_h,
    W_h = strata$W_h,
    sigma2_bh = strata$sigma2_bh,
    sigma2_wh = strata$sigma2_wh
  )
  
  return(-variance)
}