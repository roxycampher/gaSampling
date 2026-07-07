
############### TEST WITH TOY DATA ##################################

rm(list = ls())

library(dplyr)
library(readr)
library(GA)

source("variance_twostage.R")
source("cost_twostage.R")
source("check_constraints.R")
source("repair_joint.R")
source("fitness_ch_kh.R")
source("allocate_twostage_joint.R")
source("prepare_strata.R")

# Read toy data
toy_data <- read.csv("ecdc-2021-public-v2.csv")

# Convert outcome variable to numeric
toy_data$count_children_present <- readr::parse_number(
  toy_data$count_children_present
)

# Prepare strata-level data
strata <- prepare_strata(
  data = toy_data,
  strata_var = "province",
  cluster_var = "id_ward",
  unit_var = "id",
  y_var = "count_children_present",
  cost_cluster = 100,
  cost_unit = 10
)

# View prepared strata
strata

# Allocation constraints
total_n <- 900
max_cost <- 18000

# Number of strata
L <- nrow(strata)

# Bounds for c_h and k_h
lower_c <- rep(1L, L)
upper_c <- strata$C_h

lower_k <- rep(1L, L)
upper_k <- rep(total_n, L)

# ------------------------------------------------------------
# 1. Manually create an initial allocation to test repair_ckh()
# ------------------------------------------------------------

# Choose a simple starting allocation.
# This is not meant to be optimal.
# It is just used to check whether repair_ckh() works.
initial_c_h <- pmin(rep(5L, L), upper_c)
initial_k_h <- rep(10L, L)

# Check initial total sample size before repair
sum(initial_c_h * initial_k_h)

# Repair the allocation so that sum(c_h * k_h) = total_n
repaired <- repair_ckh(
  c_h = initial_c_h,
  k_h = initial_k_h,
  total_n = total_n,
  lower_c = lower_c,
  upper_c = upper_c,
  lower_k = lower_k,
  upper_k = upper_k
)

# View repaired allocation
repaired

# Check repaired total sample size
sum(repaired$c_h * repaired$k_h)

# Check repaired c_h bounds
all(repaired$c_h >= lower_c)
all(repaired$c_h <= upper_c)

# Check repaired k_h bounds
all(repaired$k_h >= lower_k)
all(repaired$k_h <= upper_k)

# Create a small allocation table for the repaired solution
repaired_allocation <- data.frame(
  stratum = strata$stratum,
  N_h = strata$N_h,
  C_h = strata$C_h,
  c_h = repaired$c_h,
  k_h = repaired$k_h,
  n_h = repaired$c_h * repaired$k_h
)

repaired_allocation

# ------------------------------------------------------------
# 2. Check the cost, variance, and constraints of repaired solution
# ------------------------------------------------------------

repaired_variance <- variance_twostage(
  c_h = repaired$c_h,
  k_h = repaired$k_h,
  W_h = strata$W_h,
  sigma2_bh = strata$sigma2_bh,
  sigma2_wh = strata$sigma2_wh
)

repaired_cost <- cost_twostage(
  c_h = repaired$c_h,
  k_h = repaired$k_h,
  cost_cluster = strata$cost_cluster,
  cost_unit = strata$cost_unit
)

repaired_check <- check_constraints(
  c_h = repaired$c_h,
  k_h = repaired$k_h,
  total_n = total_n,
  max_cost = max_cost,
  cost_cluster = strata$cost_cluster,
  cost_unit = strata$cost_unit
)

repaired_variance
repaired_cost
repaired_check

# ------------------------------------------------------------
# 3. Check fitness function using repaired allocation
# ------------------------------------------------------------

fitness_result <- fitness_ckh_constrained(
  x = c(repaired$c_h, repaired$k_h),
  strata = strata,
  total_n = total_n,
  max_cost = max_cost,
  max_c_h = strata$C_h,
  max_k_h = rep(total_n, L)
)

fitness_result

# ------------------------------------------------------------
# 4. Run the joint GA allocation
# ------------------------------------------------------------

ga_result <- allocate_twostage_ga_joint(
  strata = strata,
  total_n = total_n,
  max_cost = max_cost,
  popSize = 100,
  maxiter = 200,
  run = 50,
  seed = 123,
  suggestion = c(repaired$c_h, repaired$k_h) # start with a good allocation
)

# ------------------------------------------------------------
# 5. Inspect GA results
# ------------------------------------------------------------

ga_result$allocation
ga_result$variance
ga_result$total_n
ga_result$total_cost
ga_result$feasible
ga_result$constraint_check

# ------------------------------------------------------------
# 6. Manual checks on GA result
# ------------------------------------------------------------

# Check total sample size
sum(ga_result$allocation$c_h * ga_result$allocation$k_h)

# Check c_h does not exceed available clusters C_h
all(ga_result$allocation$c_h <= strata$C_h)

# Check c_h and k_h are at least 1
all(ga_result$allocation$c_h >= 1)
all(ga_result$allocation$k_h >= 1)

# Check cost constraint
ga_result$total_cost <= max_cost

# Check feasibility
ga_result$feasible

# ------------------------------------------------------------
# 7. Compare repaired manual allocation with GA allocation
# ------------------------------------------------------------

comparison <- data.frame(
  stratum = strata$stratum,
  C_h = strata$C_h,
  manual_c_h = repaired$c_h,
  manual_k_h = repaired$k_h,
  manual_n_h = repaired$c_h * repaired$k_h,
  ga_c_h = ga_result$allocation$c_h,
  ga_k_h = ga_result$allocation$k_h,
  ga_n_h = ga_result$allocation$c_h * ga_result$allocation$k_h
)

comparison

# Compare objective values
repaired_variance
ga_result$variance

# Compare costs
repaired_cost
ga_result$total_cost