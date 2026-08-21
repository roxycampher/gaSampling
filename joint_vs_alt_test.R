# =============================================================================
# 
# Executable script to run the optimizations and print the comparison summary
# =============================================================================
rm(list = ls())
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
source("feasible_joint_start.R")
source("allocate_twostage_joint.R")
source("alternative_joint.R")
source("frame_c.R")

# 3. Clean and process the dataset
cat("\n[1/4] Loading and cleaning ECDC dataset...")
clean_data <- read_survey_data(
  data       = read.csv("ecdc-2021-public-v2.csv"),
  stratum_id = "province",
  cluster_id = "id_ward",
  y_var      = "count_children_present",
  na_action  = "warn"
)

# 5. Define Shared Optimization Configurations
TARGET_N   <- 5000
MAX_BUDGET <- 1400000
GA_SEED    <- 42
POP_SIZE   <- 20
MAX_ITER   <- 25
RUN_LIMIT  <- 8
joint_ga_max_candidate_c <- 1000

cost_clust <- 2500
cost_unit <- 100
frame_ch <- frame_c(clean_data, MAX_BUDGET, cost_unit, cost_clust, TARGET_N)

# 4. Compute structural sampling frame rules & suggestions
cat("\n[2/4] Generating stratum metrics and initial suggestions...")
frame <- prepare_sampling_frame(
  data         = clean_data,
  y_var        = "count_children_present",
  total_c      = frame_ch,
  cost_cluster = cost_clust,  
  cost_unit    = cost_unit   
)

strata <- frame$strata

strata # check

L <- nrow(strata)

lower_c <- rep(1L, L)
upper_c <- strata$C_h

lower_k <- rep(1L, L)
upper_k <- rep(TARGET_N, L)

joint_start <- find_feasible_joint_start(
  strata = strata,
  total_n = TARGET_N,
  max_cost = MAX_BUDGET,
  lower_c = lower_c,
  upper_c = upper_c,
  lower_k = lower_k,
  upper_k = upper_k,
  max_candidate_c = min(joint_ga_max_candidate_c, frame_ch)
)

joint_starting_suggestions <- joint_start$all_starts |>
  dplyr::rename(
    repaired_sample_size_gap = sample_size_gap,
    repaired_feasible = feasible
  )

default_start <- joint_starting_suggestions |>
  dplyr::filter(starting_suggestion == "3. Neyman proportional")

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
    suggestion = c(default_start$repaired_c_h, default_start$repaired_k_h)
  )
})

# 7. Execute Proposed Decoupled Algorithm (L Dimensions)
cat("\n--- Running Alternative Decoupled GA (L Dimensions) ---\n")
time_decoupled <- system.time({
  res_decoupled <- allocate_twostage_ga_alternative(
    strata     = frame$strata,
    total_n    = TARGET_N,
    max_cost   = MAX_BUDGET,
    total_c    = 1000,           # Enforcing the global cluster cap strictly!
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
  
  cat("\n========================================================================\n")
  cat("STRATUM-LEVEL DETAILED ALLOCATION & VARIANCE COMPARISON\n")
  cat("========================================================================\n")
  cat(sprintf("%-15s | %-32s | %-32s\n", 
              "Stratum Name", 
              "Joint GA (Baseline)", 
              "Decoupled GA (Alternative)"))
  cat(sprintf("%-15s | %-32s | %-32s\n", 
              "", 
              "c_h, k_h (n_h) [B / W]", 
              "c_h, k_h (n_h) [B / W]"))
  cat("----------------|----------------------------------|----------------------------------\n")
  
  for (i in 1:nrow(res_joint$allocation)) {
    j_row <- res_joint$allocation[i, ]
    d_row <- res_decoupled$allocation[i, ]
    
    # Format strings for clean alignment
    j_str <- sprintf("c=%-2d, k=%-2d (n=%-3d) [%.3f/%.3f]", 
                     j_row$c_h, j_row$k_h, j_row$n_h, j_row$var_between, j_row$var_within)
    d_str <- sprintf("c=%-2d  , k=%-2d (n=%-3d) [%.3f/%.3f]", 
                     d_row$c_h, d_row$k_h, d_row$n_h, d_row$var_between, d_row$var_within)
    
    cat(sprintf("%-15s | %-32s | %-32s\n", j_row$stratum, j_str, d_str))
  }
  cat("========================================================================\n")
}

print_benchmark_summary(res_joint, res_decoupled, time_joint, time_decoupled)
