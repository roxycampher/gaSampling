allocate_twostage_ga_joint <- function(strata,
                                       total_n,
                                       max_cost,
                                       popSize = 100,
                                       maxiter = 1000,
                                       run = 50,
                                       seed = NULL,
                                       suggestion = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  L <- nrow(strata)
  
  if (!"C_h" %in% names(strata)) {
    stop("strata must contain C_h, the number of available clusters per stratum.")
  }
  
  lower_c <- rep(1L, L)
  lower_k <- rep(1L, L)
  upper_c <- strata$C_h
  upper_k <- rep(total_n, L)
  
  min_total_n <- sum(lower_c * lower_k)
  max_total_n <- sum(upper_c * upper_k)
  
  if (total_n < min_total_n || total_n > max_total_n) {
    stop(
      "total_n is outside the feasible range for the joint c_h and k_h bounds. ",
      "Feasible range: ", min_total_n, " to ", max_total_n, "."
    )
  }
  
  lower <- c(
    lower_c,
    lower_k
  )
  
  upper <- c(
    upper_c,
    upper_k
  )
  
  if (is.null(suggestion)) {
    
    starting_c_h <- pmax(lower_c, pmin(upper_c, round(total_n / L / 2)))
    starting_k_h <- pmax(lower_k, round(total_n / sum(starting_c_h)))
    
    starting <- repair_ckh(
      c_h = starting_c_h,
      k_h = starting_k_h,
      total_n = total_n,
      lower_c = lower_c,
      upper_c = upper_c,
      lower_k = lower_k,
      upper_k = upper_k
    )
    
    suggestions <- matrix(
      c(starting$c_h, starting$k_h),
      nrow = 1
    )
    
  } else {
    
    suggestions <- matrix(
      suggestion,
      nrow = 1
    )
  }
  
  ga_fit <- GA::ga(
    type = "real-valued",
    fitness = function(x) {
      fitness_ckh_constrained(
        x = x,
        strata = strata,
        total_n = total_n,
        max_cost = max_cost,
        min_c_h = 1L,
        min_k_h = 1L,
        max_c_h = upper_c,
        max_k_h = upper_k
      )
    },
    lower = lower,
    upper = upper,
    suggestions = suggestions,
    popSize = popSize,
    maxiter = maxiter,
    run = run
  )
  
  best_x <- ga_fit@solution[1, ]
  
  best_c_h <- round(best_x[1:L])
  best_k_h <- round(best_x[(L + 1):(2 * L)])
  
  # repair final solution
  repaired <- repair_ckh(
    c_h = best_c_h,
    k_h = best_k_h,
    total_n = total_n,
    lower_c = lower_c,
    upper_c = upper_c,
    lower_k = lower_k,
    upper_k = upper_k
  )
  
  best_c_h <- repaired$c_h
  best_k_h <- repaired$k_h

  allocation <- data.frame(
    stratum = strata$stratum,
    c_h = best_c_h,
    k_h = best_k_h,
    n_h = best_c_h * best_k_h
  ) 
  

  var_between_h <- (strata$W_h^2 * strata$sigma2_bh) / best_c_h
  var_within_h  <- (strata$W_h^2 * strata$sigma2_wh) / (best_c_h * best_k_h)
  var_total_h   <- var_between_h + var_within_h
  
  allocation$var_between <- var_between_h
  allocation$var_within  <- var_within_h
  allocation$var_total   <- var_total_h
  
  # Global sums for the summary statistics
  total_variance <- sum(var_total_h)
  total_between  <- sum(var_between_h)
  total_within   <- sum(var_within_h)
  
  # =========================================================================
  
  total_cost <- cost_twostage(
    c_h = best_c_h,
    k_h = best_k_h,
    cost_cluster = strata$cost_cluster,
    cost_unit = strata$cost_unit
  ) 
  
  constraint_check <- check_constraints(
    c_h = best_c_h,
    k_h = best_k_h,
    total_n = total_n,
    max_cost = max_cost,
    cost_cluster = strata$cost_cluster,
    cost_unit = strata$cost_unit,
    min_clusters = 1L,
    min_units = 1L,
    max_clusters = upper_c,
    max_units = upper_k
  )
  
  # --- NOTE ON FEASIBILITY GATEKEEPER ---
  # If you changed repair_ckh to handle minimums (>=), keep repaired$repaired_exactly.
  # If your repair function overshoots and returns FALSE for repaired_exactly, 
  # remove '&& repaired$repaired_exactly' so it relies strictly on constraint_check$feasible.
  is_feasible <- constraint_check$feasible && repaired$repaired_exactly
  
  list(
    feasible         = is_feasible,
    variance         = total_variance,  # Preserves keyword matching for comparison script
    between_variance = total_between,    
    within_variance  = total_within,     
    total_n          = sum(best_c_h * best_k_h),  
    total_cost       = total_cost,  
    allocation       = allocation, # Contains the new row-by-row columns
    constraint_check = constraint_check,  
    repair           = repaired,  
    ga_object        = ga_fit  
  )
}