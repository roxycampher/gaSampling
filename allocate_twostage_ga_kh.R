# Only optimizing kh and keeping ch fixed
allocate_twostage_ga_kh <- function(strata,
                                 total_n,
                                 max_cost,
                                 c_h,
                                 popSize = 100,
                                 maxiter = 1000,
                                 seed = NULL,
                                 run = 50) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  } # user can set the seed
  
  L <- nrow(strata)
  
  # bounds
  lower <- rep(1L, L)
  upper <- floor(total_n / c_h)
  
  feasibility <- check_total_n_feasible( # checks if the c_h allocation can actually produce n_total
    total_n = total_n,
    c_h = c_h,
    lower = lower,
    upper = upper
  )
  
  if (!feasibility$feasible) { # warning displayed if the c_h allocation is not feasible for n_total
    warning(feasibility$message)
  }
  
  # Simple starting solution for k_h
  starting_kh_raw <- rep(total_n / sum(c_h), L)
  
  starting_kh <- repair_kh( # repair the starting solution
    x = starting_kh_raw,
    c_h = c_h,
    total_n = total_n,
    lower = lower,
    upper = upper
  )
  
  suggestions <- matrix(starting_kh, nrow = 1) # starting solution must be in matrix form for GA
  
  ga_fit <- GA::ga(
    type = "real-valued",
    fitness = function(x) {
      fitness_kh_constrained( # fitness function (the repair is in the fitness function as well)
        x = x,
        strata = strata,
        total_n = total_n,
        max_cost = max_cost,
        c_h = c_h
      )
    },
    lower = lower, # min value allowed for each gene
    upper = upper, # max value allowed for each gene
    suggestions = suggestions,
    popSize = popSize, # number of candidate solutions per generation
    maxiter = maxiter, # number of generations
    run = run # early stopping rule
  )
  
  best_kh <- repair_kh( # post repair in case the GA solution is still real valued
    x = ga_fit@solution[1, ],
    c_h = c_h,
    total_n = total_n,
    lower = lower,
    upper = upper
  )
  
  best_nh <- best_kh * c_h # final stratum sample sizes
  
  constraint_check <- check_constraints( # verify the final solution
    c_h = c_h,
    k_h = best_kh,
    total_n = total_n,
    max_cost = max_cost,
    cost_cluster = strata$cost_cluster,
    cost_unit = strata$cost_unit
  )
  
  if (!constraint_check$feasible) { # stop if no feasible allocation is found
    stop("The GA did not find a feasible allocation. Try increasing popSize/maxiter or relaxing the exact total_n constraint.")
  }
  
  # final variance calculation
  best_variance <- variance_twostage(
    c_h = c_h,
    k_h = best_kh,
    W_h = strata$W_h,
    sigma2_bh = strata$sigma2_bh,
    sigma2_wh = strata$sigma2_wh
  )
  
  # final cost calculation
  best_cost <- cost_twostage(
    c_h = c_h,
    k_h = best_kh,
    cost_cluster = strata$cost_cluster,
    cost_unit = strata$cost_unit
  )
  
  # The Output table
  allocation <- data.frame(
    stratum = strata$stratum,
    c_h = c_h,
    k_h = best_kh,
    n_h = best_nh
  )
  
  # return results
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