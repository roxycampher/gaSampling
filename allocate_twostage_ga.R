allocate_twostage_ga <- function(strata,
                                 total_n,
                                 max_cost,
                                 c_h,
                                 popSize = 100,
                                 maxiter = 1000,
                                 seed = NULL) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  L <- nrow(strata)
  
  lower <- rep(1L, L)
  upper <- floor(total_n / c_h)
  
  feasibility <- check_total_n_feasible(
    total_n = total_n,
    c_h = c_h,
    lower = lower,
    upper = upper
  )
  
  if (!feasibility$feasible) {
    warning(feasibility$message)
  }
  
  # Simple starting solution
  starting_kh_raw <- rep(total_n / sum(c_h), L)
  
  starting_kh <- repair_kh(
    x = starting_kh_raw,
    c_h = c_h,
    total_n = total_n,
    lower = lower,
    upper = upper
  )
  
  suggestions <- matrix(starting_kh, nrow = 1)
  
  ga_fit <- GA::ga(
    type = "real-valued",
    fitness = function(x) {
      fitness_kh_constrained(
        x = x,
        strata = strata,
        total_n = total_n,
        max_cost = max_cost,
        c_h = c_h
      )
    },
    lower = lower,
    upper = upper,
    suggestions = suggestions,
    popSize = popSize,
    maxiter = maxiter,
    run = 50
  )
  
  best_kh <- repair_kh(
    x = ga_fit@solution[1, ],
    c_h = c_h,
    total_n = total_n,
    lower = rep(1L, length(c_h)),
    upper = floor(total_n / c_h)
  )
  
  best_nh <- best_kh * c_h
  
  constraint_check <- check_constraints(
    c_h = c_h,
    k_h = best_kh,
    total_n = total_n,
    max_cost = max_cost,
    cost_cluster = strata$cost_cluster,
    cost_unit = strata$cost_unit
  )
  
  if (!constraint_check$feasible) {
    stop("The GA did not find a feasible allocation. Try increasing popSize/maxiter or relaxing the exact total_n constraint.")
  }
  
  best_variance <- variance_twostage(
    c_h = c_h,
    k_h = best_kh,
    W_h = strata$W_h,
    sigma2_bh = strata$sigma2_bh,
    sigma2_wh = strata$sigma2_wh
  )
  
  best_cost <- cost_twostage(
    c_h = c_h,
    k_h = best_kh,
    cost_cluster = strata$cost_cluster,
    cost_unit = strata$cost_unit
  )
  
  allocation <- data.frame(
    stratum = strata$stratum,
    c_h = c_h,
    k_h = best_kh,
    n_h = best_nh
  )
  
  list(
    allocation = allocation,
    variance = best_variance,
    total_n = sum(best_nh),
    total_cost = best_cost,
    feasible = constraint_check$feasible,
    constraint_check = constraint_check,
    ga_object = ga_fit
  )
}