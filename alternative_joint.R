# =============================================================================
# gasampling: Decoupled Two-Stage Optimal Allocation via GA
# File Name: alternative_joint.R
# =============================================================================

allocate_twostage_ga_alternative <- function(strata,
                                             total_n,
                                             max_cost,
                                             total_c,       
                                             popSize = 50,
                                             maxiter = 100,
                                             run = 20,
                                             seed = NULL,
                                             suggestion = NULL) {
  L <- nrow(strata)
  min_c <- rep(1L, L)
  max_c <- strata$C_h
  
  # =========================================================================
  # STEP 1: Compute TRUE Analytical k_h and adjust for the sample size target
  # =========================================================================
  opt_k <- numeric(L)
  for (h in 1:L) {
    val <- sqrt((strata$sigma2_wh[h] * strata$cost_cluster[h]) / 
                  (strata$sigma2_bh[h] * strata$cost_unit[h]))
    opt_k[h] <- round(val)
  }
  
  # If the total cluster limit multiplied by the optimal k_h cannot reach total_n,
  # scale k_h up uniformly up front so that achieving total_n is mathematically possible.
  max_possible_n <- total_c * mean(opt_k)
  if (max_possible_n < total_n) {
    scaling_factor <- total_n / max_possible_n
    opt_k <- round(opt_k * scaling_factor)
  }
  
  # Ensure k_h values stay within absolute boundaries
  opt_k <- pmax(opt_k, 1L)
  opt_k <- pmin(opt_k, total_n)
  
  # --- DECOUPLED FITNESS FUNCTION (L-DIMENSIONS) ---
  fitness_decoupled <- function(x) {
    c_h <- round(x)
    c_h <- pmax(c_h, min_c)
    c_h <- pmin(c_h, max_c)
    
    # Strict Repair for Global Cluster Limit: sum(c_h) <= total_c
    if (sum(c_h) > total_c) {
      while (sum(c_h) > total_c) {
        reducible_indices <- which(c_h > min_c)
        if (length(reducible_indices) == 0) break
        idx <- reducible_indices[which.max(c_h[reducible_indices])]
        c_h[idx] <- c_h[idx] - 1L
      }
    }
    
    # Evaluate using the target-aware optimal k_h values
    check <- check_constraints(
      c_h          = c_h,
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
    
    if (!check$feasible) {
      return(-1e10)
    }
    
    variance <- variance_twostage(c_h, opt_k, strata$W_h, strata$sigma2_bh, strata$sigma2_wh)
    return(-variance)
  }
  
  # --- RUN GENETIC ALGORITHM ---
  if (!is.null(seed)) set.seed(seed)
  
  ga_res <- ga(
    type        = "real-valued",
    fitness     = fitness_decoupled,
    lower       = min_c,
    upper       = max_c,
    popSize     = popSize,
    maxiter     = maxiter,
    run         = run,
    suggestions = if (!is.null(suggestion)) matrix(suggestion, nrow = 1) else NULL,
    monitor     = FALSE
  )
  
  # --- EXTRACT BEST SOLUTION ---
  best_c <- round(ga_res@solution[1, ])
  
  if (sum(best_c) > total_c) {
    while (sum(best_c) > total_c) {
      reducible_indices <- which(best_c > min_c)
      if (length(reducible_indices) == 0) break
      idx <- reducible_indices[which.max(best_c[reducible_indices])]
      best_c[idx] <- best_c[idx] - 1L
    }
  }
  
  final_check <- check_constraints(
    c_h = best_c, k_h = opt_k, total_n = total_n, max_cost = max_cost,
    cost_cluster = strata$cost_cluster, cost_unit = strata$cost_unit,
    min_clusters = min_c, min_units = rep(1L, L), max_clusters = max_c, max_units = rep(total_n, L)
  )
  
  # --- EXTRACT BEST SOLUTION & RE-RUN CONSTRAINTS FOR OUTPUT ---
  # ... (Keep your existing best_c and best_k extraction code here) ... 
  
  # Calculate stratum-level variance components explicitly
  var_between_h <- (strata$W_h^2 * strata$sigma2_bh) / best_c
  var_within_h  <- (strata$W_h^2 * strata$sigma2_wh) / (best_c * best_k)
  var_total_h   <- var_between_h + var_within_h
  
  # Sum across all strata to get global totals
  total_variance <- sum(var_total_h)
  total_between  <- sum(var_between_h)
  total_within   <- sum(var_within_h)
  
  # Return matching structured list now including variance components
  list(
    feasible        = final_check$feasible, 
    variance        = total_variance,       # Keeps compatibility with old print functions 
    between_variance = total_between,        # New metric!
    within_variance  = total_within,         # New metric!
    total_cost      = final_check$total_cost, 
    total_n         = final_check$total_n, 
    total_c         = sum(best_c), 
    allocation      = data.frame(
      stratum        = strata$stratum, 
      c_h            = best_c, 
      k_h            = best_k, 
      n_h            = final_check$n_h,
      var_between    = var_between_h,       # Stratum level breakdown
      var_within     = var_within_h,        # Stratum level breakdown
      var_total      = var_total_h
    )
  )

}