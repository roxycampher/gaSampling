# =============================================================================
# 
# Executable script to run the optimizations and print the comparison summary
# =============================================================================

# 1. Load Required External Packages
library(GA)
library(dplyr)

# 2. Source all the underlying infrastructure files
source("read_survey_data.R")
source("sampling_frame.R")
source("repair_joint.R")
source("check_constraints.R")
source("variance_twostage.R")
source("cost_twostage.R")
source("fitness_ch_kh.R")
source("allocate_twostage_joint.R")
source("alternative_joint.R")

# 3. Clean and process the dataset
cat("\n[1/4] Loading and cleaning ECDC dataset...")
clean_data <- read_survey_data(
  data       = read.csv("ecdc-2021-public-v2.csv"),
  stratum_id = "province",
  cluster_id = "dc",
  y_var      = "count_register_all",
  na_action  = "warn"
)

# 4. Compute structural sampling frame rules & suggestions
cat("\n[2/4] Generating stratum metrics and initial suggestions...")
frame <- prepare_sampling_frame(
  data         = clean_data,
  y_var        = "count_register_all",
  total_c      = 30,
  cost_cluster = 500,  
  cost_unit    = 50    
)

# 5. Define Shared Optimization Configurations
TARGET_N   <- 1000
MAX_BUDGET <- 70000
GA_SEED    <- 42
POP_SIZE   <- 80
MAX_ITER   <- 300
RUN_LIMIT  <- 30

# 6. Execute Baseline Joint Algorithm (2L Dimensions)
cat("\n[3/4] Running Baseline Joint GA (2L Dimensions)... Please wait.\n")
time_joint <- system.time({
  res_joint <- allocate_twostage_ga_joint(
    strata     = frame$strata,
    total_n    = TARGET_N,
    max_cost   = MAX_BUDGET,
    popSize    = POP_SIZE,
    maxiter    = MAX_ITER,
    run        = RUN_LIMIT,
    seed       = GA_SEED,
    suggestion = c(frame$c_h, rep(round(TARGET_N / sum(frame$c_h)), nrow(frame$strata))) 
  )
})

# 7. Execute Proposed Decoupled Algorithm (L Dimensions)
cat("\n--- Running Alternative Decoupled GA (L Dimensions) ---\n")
time_decoupled <- system.time({
  res_decoupled <- allocate_twostage_ga_alternative(
    strata     = frame$strata,
    total_n    = TARGET_N,
    max_cost   = MAX_BUDGET,
    total_c    = 30,           # Enforcing the global cluster cap strictly!
    popSize    = POP_SIZE,
    maxiter    = MAX_ITER,
    run        = RUN_LIMIT,
    seed       = GA_SEED,
    suggestion = frame$c_h  
  )
})


print_benchmark_summary <- function(res_joint, res_decoupled, time_joint, time_decoupled) {
  cat("\n========================================================================\n")
  cat("TWO-STAGE SAMPLING GENETIC ALGORITHM OPTIMIZATION SUMMARY\n")
  cat("========================================================================\n")
  
  cat(sprintf("\nBaseline Model: Joint GA (2L Dimension Search Space)\n"))
  cat(sprintf("  Feasibility Status     : %s\n", if(res_joint$feasible) "FEASIBLE" else "INFEASIBLE"))
  cat(sprintf("  Total Wall-Clock Time  : %.4f seconds\n", time_joint["elapsed"]))
  cat(sprintf("  Final Total Cost       : R %.2f\n", res_joint$total_cost))
  cat(sprintf("  Achieved Sample Size   : %d\n", res_joint$total_n))
  cat(sprintf("  Total Variance (V)     : %.8f\n", res_joint$variance))
  cat(sprintf("    └─ Between-Strata Var: %.8f\n", res_joint$between_variance))
  cat(sprintf("    └─ Within-Strata Var : %.8f\n", res_joint$within_variance))
  
  cat(sprintf("\nProposed Alternative Model: Decoupled GA (L Dimension Search Space)\n"))
  cat(sprintf("  Feasibility Status     : %s\n", if(res_decoupled$feasible) "FEASIBLE" else "INFEASIBLE"))
  cat(sprintf("  Total Wall-Clock Time  : %.4f seconds\n", time_decoupled["elapsed"]))
  cat(sprintf("  Final Total Cost       : R %.2f\n", res_decoupled$total_cost))
  cat(sprintf("  Achieved Sample Size   : %d\n", res_decoupled$total_n))
  cat(sprintf("  Total Variance (V)     : %.8f\n", res_decoupled$variance))
  cat(sprintf("    └─ Between-Strata Var: %.8f\n", res_decoupled$between_variance))
  cat(sprintf("    └─ Within-Strata Var : %.8f\n", res_decoupled$within_variance))
  
  # Structural math metrics
  time_pct <- ((time_joint["elapsed"] - time_decoupled["elapsed"]) / time_joint["elapsed"]) * 100
  var_pct  <- ((res_joint$variance - res_decoupled$variance) / res_joint$variance) * 100
  
  cat("\nAlgorithm Efficiency Comparison:\n")
  cat(sprintf("  Wall-Clock Time Reduction           : +%.2f%%\n", time_pct))
  cat(sprintf("  Target Variance Minimization Delta  : %+.2f%% (Alternative is %s)\n", 
              var_pct, if(var_pct > 0) "better" else "worse"))
  cat("========================================================================\n")
  cat("\n========================================================================\n")
  cat("STRATUM-LEVEL VARIANCE COMPONENT BREAKDOWN COMPARISON\n")
  cat("========================================================================\n")
  cat(sprintf("%-16s | %-30s | %-30s\n", "Stratum Name", "Joint GA (Between / Within)", "Decoupled GA (Between / Within)"))
  cat("-----------------|--------------------------------|--------------------------------\n")
  
  for (i in 1:nrow(frame$strata)) {
    joint_row     <- res_joint$allocation[i, ]
    decoupled_row <- res_decoupled$allocation[i, ]
    
    cat(sprintf("%-16s | B: %.4f / W: %.4f (n=%-4d) | B: %.4f / W: %.4f (n=%-4d)\n",
                joint_row$stratum,
                joint_row$var_between, joint_row$var_within, joint_row$n_h,
                decoupled_row$var_between, decoupled_row$var_within, decoupled_row$n_h))
  }
  cat("========================================================================\n")
}

print_benchmark_summary(res_joint, res_decoupled, time_joint, time_decoupled)
