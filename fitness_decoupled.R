fitness_decoupled_clean <- function(x,
                                    opt_k,
                                    strata,
                                    total_n,
                                    max_cost,
                                    min_c_h = 1L,
                                    penalty_weight = 1e5) {
  
  L <- nrow(strata)
  
  # Decode c_h from decision vector x (L dimensions)
  c_h <- round(x)
  c_h <- pmax(c_h, min_c_h)
  c_h <- pmin(c_h, strata$C_h)
  
  # 1. Compute Total Cost and Achieved Sample Size directly
  current_cost <- sum(strata$cost_cluster * c_h + strata$cost_unit * c_h * opt_k)
  achieved_n   <- sum(c_h * opt_k)
  
  # 2. Compute Variance
  variance <- variance_twostage(
    c_h       = c_h,
    k_h       = opt_k,
    W_h       = strata$W_h,
    sigma2_bh = strata$sigma2_bh,
    sigma2_wh = strata$sigma2_wh
  )
  
  # 3. Calculate Smooth Penalties (No early hard return)
  cost_overrun     <- max(0, current_cost - max_cost)
  sample_shortfall <- max(0, total_n - achieved_n)
  
  # Normalize penalties so they scale predictably with variance
  penalty <- penalty_weight * (cost_overrun / max_cost) + 
    penalty_weight * (sample_shortfall / total_n)
  
  # GA maximizes fitness, so return negative (variance + penalty)
  return(-(variance + penalty))
}