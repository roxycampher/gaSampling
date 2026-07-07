# Optimizing both kh and ch
allocate_twostage_ga_joint <- function(strata,
                                       total_n,
                                       max_cost,
                                       popSize = 100,
                                       maxiter = 1000,
                                       run = 50,
                                       seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  L <- nrow(strata)
  
  lower_c <- rep(1L, L)
  lower_k <- rep(1L, L)
  upper_c <- strata$M_h
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
    rep(1L, L),  # lower bounds for c_h
    rep(1L, L)   # lower bounds for k_h
  )
  
  upper <- c(
    strata$M_h,                  # max possible clusters in each stratum
    rep(total_n, L)              # rough upper bound for k_h
  )
  
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
    constraint_check = constraint_check,
    repair = repaired,
    ga_fit = ga_fit
  )
}