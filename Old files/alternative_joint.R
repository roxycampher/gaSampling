# =============================================================================
# gasampling: Decoupled Two-Stage Optimal Allocation via GA
# File Name: alternative_joint.R
# =============================================================================

allocate_twostage_ga_alternative <- function(strata,
                                             total_n,
                                             max_cost,
                                             popSize = 100,      # Increased population size for better search
                                             maxiter = 500,      # High maxiter
                                             run = 50,           # Higher stall limit before stopping
                                             seed = NULL) {      # Removed suggestion parameter
  L <- nrow(strata)
  min_c <- rep(1L, L)
  max_c <- strata$C_h
  
  # =========================================================================
  # STEP 1: Compute Constrained Analytical k_h
  # =========================================================================
  opt_k <- numeric(L)
  for (h in 1:L) {
    val <- sqrt((strata$sigma2_wh[h] * strata$cost_cluster[h]) / 
                  (strata$sigma2_bh[h] * strata$cost_unit[h]))
    opt_k[h] <- val
  }
  
  # Ensure k_h is scaled sufficiently high to realistically achieve target_n
  cost_per_cluster <- strata$cost_cluster + strata$cost_unit * opt_k
  est_total_c <- max_cost / mean(cost_per_cluster)
  est_raw_n <- est_total_c * mean(opt_k)
  
  if (est_raw_n < total_n) {
    scale_factor <- total_n / est_raw_n
    opt_k <- opt_k * scale_factor
  }
  
  opt_k <- round(opt_k)
  opt_k <- pmax(opt_k, 1L)
  
  # =========================================================================
  # STEP 2: Dedicated Decoupled Fitness Function
  # =========================================================================
  fitness_decoupled_clean <- function(x) {
    c_h <- round(x)
    c_h <- pmax(c_h, min_c)
    c_h <- pmin(c_h, max_c)
    
    current_cost <- sum(strata$cost_cluster * c_h + strata$cost_unit * c_h * opt_k)
    achieved_n   <- sum(c_h * opt_k)
    
    variance <- variance_twostage(
      c_h       = c_h,
      k_h       = opt_k,
      W_h       = strata$W_h,
      sigma2_bh = strata$sigma2_bh,
      sigma2_wh = strata$sigma2_wh
    )
    
    # Quadratic relative penalties for smooth gradient signals
    cost_overrun     <- max(0, (current_cost - max_cost) / max_cost)^2
    sample_shortfall <- max(0, (total_n - achieved_n) / total_n)^2
    
    # Penalty multiplier scaled to variance magnitude (~1.0)
    penalty <- 100 * cost_overrun + 100 * sample_shortfall
    
    return(-(variance + penalty))
  }
  
  # =========================================================================
  # STEP 3: Run GA Search
  # =========================================================================
  if (!is.null(seed)) set.seed(seed)
  
  ga_res <- ga(
    type        = "real-valued",
    fitness     = fitness_decoupled_clean,
    lower       = min_c,
    upper       = max_c,
    popSize     = popSize,
    maxiter     = maxiter,
    run         = run,
    monitor     = FALSE,
    seed        = seed
  )
  
  # =========================================================================
  # STEP 4: Extract Best Solution
  # =========================================================================
  best_c <- round(ga_res@solution[1, ])
  best_c <- pmax(best_c, min_c)
  best_c <- pmin(best_c, max_c)
  
  final_check <- check_constraints(
    c_h          = best_c, 
    k_h          = opt_k, 
    total_n      = total_n, 
    max_cost     = max_cost,
    cost_cluster = strata$cost_cluster, 
    cost_unit    = strata$cost_unit,
    min_clusters = min_c, 
    min_units    = rep(1L, L), 
    max_clusters = max_c, 
    max_units    = rep(total_n, L)
  )
  
  var_between_h <- (strata$W_h^2 * strata$sigma2_bh) / best_c
  var_within_h  <- (strata$W_h^2 * strata$sigma2_wh) / (best_c * opt_k)
  var_total_h   <- var_between_h + var_within_h
  
  list(
    feasible         = final_check$feasible, 
    variance         = sum(var_total_h),       
    between_variance = sum(var_between_h),        
    within_variance  = sum(var_within_h),         
    total_cost       = final_check$total_cost, 
    total_n          = final_check$total_n, 
    total_c          = sum(best_c), 
    allocation       = data.frame(
      stratum        = strata$stratum, 
      c_h            = best_c, 
      k_h            = opt_k, 
      n_h            = final_check$n_h,
      var_between    = var_between_h,       
      var_within     = var_within_h,        
      var_total      = var_total_h
    )
  )
}