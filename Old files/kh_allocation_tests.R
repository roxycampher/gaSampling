rm(list = ls())
library(dplyr)
library(readr)
library(GA)

source("variance_twostage.R")
source("cost_twostage.R")
source("check_constraints.R")
source("repair_kh.R")
source("sample_size_feasibility.R")
source("fitness_kh.R")
source("prepare_strata.R")
source("allocate_twostage_ga_kh.R")

# The data cannot be stored in the github as it is to large - will just remove before commits and pushing
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

(total_n <- sum(c_h * k_h))
max_cost <- 20000

# Test variance function
var_result <- variance_twostage(
  c_h = c_h,
  k_h = k_h,
  W_h = strata$W_h,
  sigma2_bh = strata$sigma2_bh,
  sigma2_wh = strata$sigma2_wh
)

var_result

# Test cost function
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

# Check fitness equals negative variance (We use - variance because we will use GA to maximize the -variance)
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
source("fitness_kh.R")
source("prepare_strata.R")
source("allocate_twostage_ga_kh.R")

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

ga_result <- allocate_twostage_ga_kh(
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
library(ggplot2)
seeds <- 1:25

# run the GA over seeds
ga_runs <- lapply(seeds, function(s) {
  allocate_twostage_ga_kh(
    strata = strata,
    total_n = 900,
    max_cost = 18000,
    c_h = c_h,
    popSize = 100,
    maxiter = 1000,
    seed = s
  )
})

# Take the list of GA outputs and then extract the summary statistics we need
run_summary <- data.frame(
  seed = seeds,
  variance = sapply(ga_runs, function(x) x$variance),
  total_n = sapply(ga_runs, function(x) x$total_n),
  total_cost = sapply(ga_runs, function(x) x$total_cost),
  feasible = sapply(ga_runs, function(x) x$feasible)
)

run_summary

min(run_summary$variance) # best variance achieved
mean(run_summary$variance)
sd(run_summary$variance) # should be small (stable)
sum(run_summary$feasible) # how many feasible solutions are there?

run_summary %>%
  ggplot( aes(x=seeds, y=variance)) +
  geom_line() +
  geom_point() +
  ggtitle("Variance obtained across different seeds.")

# Best GA run
best_index <- which.min(run_summary$variance)

best_run <- ga_runs[[best_index]] # GA output for the best run

best_run$allocation
best_run$variance
best_run$total_n
best_run$total_cost
best_run$feasible

# Average allocations across runs (how stable allocations are across stratum)
allocations <- dplyr::bind_rows(
  lapply(seq_along(ga_runs), function(i) {
    ga_runs[[i]]$allocation |>
      dplyr::mutate(seed = seeds[i])
  })
)

allocation_summary <- allocations |>
  dplyr::group_by(stratum) |>
  dplyr::summarise(
    mean_k_h = mean(k_h),
    sd_k_h = sd(k_h),
    min_k_h = min(k_h), # smallest k_h across all runs
    max_k_h = max(k_h), # largest k_h across all runs
    mean_n_h = mean(n_h), # avg stratum sample size across the GA runs
    .groups = "drop"
  )

allocation_summary

#### Comparison to equal allocation ###

# EQUAL ALLOCATION
k_equal <- rep(10, nrow(strata)) # set an equal allocation for k_h across strata

variance_equal <- variance_twostage(
  c_h = c_h, # fixed the same as before
  k_h = k_equal, # equal k_h allocation
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

# COMPARISON 
best_run$variance
variance_equal

better1 <- best_run$variance < variance_equal
better1

### Compare with proportional allocation ###

# Raw proportional stratum sample sizes
n_prop_raw <- total_n * strata$W_h

k_prop_raw <- n_prop_raw / c_h

# make integer
k_prop <- repair_kh(
  x = k_prop_raw,
  c_h = c_h,
  total_n = total_n,
  lower = rep(1L, length(c_h)),
  upper = floor(total_n / c_h)
)

n_prop <- c_h * k_prop

# calculate variance
variance_prop <- variance_twostage(
  c_h = c_h,
  k_h = k_prop,
  W_h = strata$W_h,
  sigma2_bh = strata$sigma2_bh,
  sigma2_wh = strata$sigma2_wh
)

# calculate cost
cost_prop <- cost_twostage(
  c_h = c_h,
  k_h = k_prop,
  cost_cluster = strata$cost_cluster,
  cost_unit = strata$cost_unit
)

# check feasibility
check_prop <- check_constraints(
  k_h = k_prop,
  c_h = c_h,
  total_n = total_n,
  max_cost = max_cost,
  cost_cluster = strata$cost_cluster,
  cost_unit = strata$cost_unit
)

# create proportional allocation table
allocation_prop <- data.frame(
  stratum = strata$stratum,
  W_h = strata$W_h,
  c_h = c_h,
  k_h = k_prop,
  n_h = n_prop,
  n_prop_raw = n_prop_raw
)

allocation_prop
variance_prop
cost_prop
check_prop

# compare
variance_prop
best_run$variance

better2 <- best_run$variance < variance_prop
better2

########## TEST WHEN WE HAVE PROP ALLOC FOR CH ################################
# Raw proportional cluster allocation
total_c <- 90

c_h_raw <- total_c * strata$W_h

# Round to integers
c_h_prop <- round(c_h_raw)

ga_result_prop_ch <- allocate_twostage_ga_kh(
  strata = strata,
  total_n = total_n,
  max_cost = max_cost,
  c_h = c_h_prop,
  popSize = 100,
  maxiter = 500,
  seed = 123
)

ga_result_prop_ch$allocation
ga_result_prop_ch$variance
ga_result_prop_ch$total_n
ga_result_prop_ch$total_cost
ga_result_prop_ch$feasible
ga_result_prop_ch$constraint_check

