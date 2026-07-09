rm(list = ls())

source("read_survey_data.R")
source("sampling_frame.R")
source("variance_twostage.R")
source("cost_twostage.R")
source("check_constraints.R")
source("sample_size_feasibility.R")
source("repair_kh.R")
source("fitness_kh.R")
source("allocate_twostage_ga_kh.R")

ecdc_ssu <- read.csv("ecdc-2021-public-v2.csv")

clean_data <- read_survey_data(
  data       = ecdc_ssu,
  stratum_id = "province",
  cluster_id = "id_ward",
  y_var      = "count_children_present",
  na_action  = "warn"
)

# Prepare frame from SSU-level data
frame <- prepare_sampling_frame(
  data         = clean_data,
  y_var        = "count_children_present",
  total_c      = 90,
  cost_cluster = 100,
  cost_unit    = 10
)

# Extract GA inputs
strata <- frame$strata
c_h <- rep(10, nrow(strata))

# Optimise k_h while holding c_h fixed
kh_result <- allocate_twostage_ga_kh(
  strata = strata,
  c_h = c_h,
  total_n = 900,
  max_cost = 20000,
  popSize = 100,
  maxiter = 200,
  seed = 123
)

kh_result$allocation
kh_result$variance
kh_result$total_cost
kh_result$feasible

# Final allocation sample sizes
final_allocation <- kh_result$allocation

final_allocation$n_h <- final_allocation$c_h * final_allocation$k_h

print(final_allocation)

