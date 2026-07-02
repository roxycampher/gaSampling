# Optimizing both kh and ch
allocate_twostage_ga_joint <- function(strata,
                                       total_n,
                                       max_cost,
                                       popSize = 100,
                                       maxiter = 1000,
                                       seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  L <- nrow(strata)
  
  lower <- c(
    rep(1L, L),  # lower bounds for c_h
    rep(1L, L)   # lower bounds for k_h
  )
  
  upper <- c(
    strata$M_h,                  # max possible clusters in each stratum
    rep(total_n, L)              # rough upper bound for k_h
  )
  
  ga_fit <- GA::ga(
    type = "real-valued",
    fitness = function(x) {
      fitness_ckh_constrained(
        x = x,
        strata = strata,
        total_n = total_n,
        max_cost = max_cost
      )
    },
    lower = lower,
    upper = upper,
    popSize = popSize,
    maxiter = maxiter,
    run = 50
  )
  
  best_x <- ga_fit@solution[1, ]
  
  best_c_h <- round(best_x[1:L])
  best_k_h <- round(best_x[(L + 1):(2 * L)])
  
  # repair final solution (need new repair function for joint)
  repaired <- repair_ckh(
    c_h = best_c_h,
    k_h = best_k_h,
    total_n = total_n,
    lower_c = rep(1L, L),
    lower_k = rep(1L, L),
    upper_c = strata$M_h,
    upper_k = rep(total_n, L)
  ) # proposed structure
  
  best_c_h <- repaired$c_h
  best_k_h <- repaired$k_h
  
  constraint_check <- check_constraints(
    c_h = best_c_h,
    k_h = best_k_h,
    total_n = total_n,
    max_cost = max_cost,
    cost_cluster = strata$cost_cluster,
    cost_unit = strata$cost_unit
  )
  
  allocation <- data.frame(
    stratum = strata$stratum,
    c_h = best_c_h,
    k_h = best_k_h,
    n_h = best_c_h * best_k_h
  )
  
  variance <- variance_twostage(
    c_h = best_c_h,
    k_h = best_k_h,
    W_h = strata$W_h,
    sigma2_bh = strata$sigma2_bh,
    sigma2_wh = strata$sigma2_wh
  )
  
  total_cost <- cost_twostage(
    c_h = best_c_h,
    k_h = best_k_h,
    cost_cluster = strata$cost_cluster,
    cost_unit = strata$cost_unit
  )
  
  list(
    allocation = allocation,
    variance = variance,
    total_n = sum(best_c_h * best_k_h),
    total_cost = total_cost,
    feasible = constraint_check$feasible,
    ga_fit = ga_fit
  )
}