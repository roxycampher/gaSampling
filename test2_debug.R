library(dplyr)
library(readr)
library(GA)

source("variance_twostage.R")
source("cost_twostage.R")
source("check_constraints.R")
source("repair_kh.R")
source("sample_size_feasibility.R")
source("fitness.R")
source("prepare_strata.R")
source("allocate_twostage_ga.R")

toy_data <- read.csv("ecdc-2021-public-v2.csv")

# Clean outcome variable
toy_data$count_children_present <- readr::parse_number(
  toy_data$count_children_present
)

# Prepare strata from toy_data
strata <- prepare_strata(
  data = toy_data,
  strata_var = "province",
  cluster_var = "id_ward",
  unit_var = "id",
  y_var = "count_children_present",
  cost_cluster = 100,
  cost_unit = 10
)

# Basic checks
strata
str(strata)
colSums(is.na(strata))
sum(strata$W_h)

# Candidate allocation
c_h <- rep(10, nrow(strata))
k_h <- rep(10, nrow(strata)) # we fit both to test the functions - not an accurate final result just a test

total_n <- sum(c_h * k_h)
max_cost <- 20000

# Test variance
var_result <- variance_twostage(
  c_h = c_h,
  k_h = k_h,
  W_h = strata$W_h,
  sigma2_bh = strata$sigma2_bh,
  sigma2_wh = strata$sigma2_wh
)

var_result

# Test cost
cost_result <- cost_twostage(
  c_h = c_h,
  k_h = k_h,
  cost_cluster = strata$cost_cluster,
  cost_unit = strata$cost_unit
)

cost_result

# Test constraints
constraint_result <- check_constraints(
  c_h = c_h,
  k_h = k_h,
  total_n = total_n,
  max_cost = max_cost,
  cost_cluster = strata$cost_cluster,
  cost_unit = strata$cost_unit
)

constraint_result

# Test fitness
fitness_result <- fitness_kh_constrained(
  x = k_h,
  strata = strata,
  total_n = total_n,
  max_cost = max_cost,
  c_h = c_h
)

fitness_result

# Check fitness equals negative variance (We use - variance because we will use GA to maximise the -variance)
all.equal(fitness_result, -var_result)


########################### We will now test the allocation on the toy_data ############################
rm(list = ls())

library(dplyr)
library(readr)
library(GA)

source("variance_twostage.R")
source("cost_twostage.R")
source("check_constraints.R")
source("repair_kh.R")
source("sample_size_feasibility.R")
source("fitness.R")
source("prepare_strata.R")
source("allocate_twostage_ga.R")

toy_data <- read.csv("ecdc-2021-public-v2.csv")

toy_data$count_children_present <- readr::parse_number(
  toy_data$count_children_present
)

strata <- prepare_strata(
  data = toy_data,
  strata_var = "province",
  cluster_var = "id_ward",
  unit_var = "id",
  y_var = "count_children_present",
  cost_cluster = 100,
  cost_unit = 10
)

c_h <- rep(10, nrow(strata)) # fix c_h

total_n <- 900
max_cost <- 20000

ga_result <- allocate_twostage_ga(
  strata = strata,
  total_n = total_n,
  max_cost = max_cost,
  c_h = c_h,
  popSize = 100,
  maxiter = 200,
  seed = 123
)

ga_result$allocation
ga_result$variance
ga_result$total_n
ga_result$total_cost
ga_result$feasible


################# TEST MULTIPLE RUNS ###########################################

set.seed(2026)

seeds <- 1:50

ga_runs <- lapply(seeds, function(s) {
  allocate_twostage_ga(
    strata = strata,
    total_n = 900,
    max_cost = 18000,
    c_h = rep(10, nrow(strata)),
    popSize = 100,
    maxiter = 1000,
    seed = s
  )
})

run_summary <- data.frame(
  seed = seeds,
  variance = sapply(ga_runs, function(x) x$variance),
  total_n = sapply(ga_runs, function(x) x$total_n),
  total_cost = sapply(ga_runs, function(x) x$total_cost),
  feasible = sapply(ga_runs, function(x) x$feasible)
)

run_summary

summary(run_summary$variance)

min(run_summary$variance)
mean(run_summary$variance)
sd(run_summary$variance) # should be small
sum(run_summary$feasible)

# Best GA run
best_index <- which.min(run_summary$variance)

best_run <- ga_runs[[best_index]]

best_run$allocation
best_run$variance
best_run$total_n
best_run$total_cost
best_run$feasible

# Average allocations across runs
allocations <- do.call(rbind, lapply(seq_along(ga_runs), function(i) {
  out <- ga_runs[[i]]$allocation
  out$seed <- seeds[i]
  out
}))

allocation_summary <- allocations |>
  dplyr::group_by(stratum) |>
  dplyr::summarise(
    mean_k_h = mean(k_h),
    sd_k_h = sd(k_h),
    min_k_h = min(k_h),
    max_k_h = max(k_h),
    mean_n_h = mean(n_h),
    .groups = "drop"
  )

allocation_summary

# Comparison to equal allocation
k_equal <- rep(10, nrow(strata))

c_h <- rep(10, nrow(strata))

variance_equal <- variance_twostage(
  c_h = c_h,
  k_h = k_equal,
  W_h = strata$W_h,
  sigma2_bh = strata$sigma2_bh,
  sigma2_wh = strata$sigma2_wh
)

cost_equal <- cost_twostage(
  c_h = c_h,
  k_h = k_equal,
  cost_cluster = strata$cost_cluster,
  cost_unit = strata$cost_unit
)

check_equal <- check_constraints(
  k_h = k_equal,
  c_h = c_h,
  total_n = 900,
  max_cost = 18000,
  cost_cluster = strata$cost_cluster,
  cost_unit = strata$cost_unit
)

variance_equal
cost_equal
check_equal

# Comnpare 
best_run$variance
variance_equal
