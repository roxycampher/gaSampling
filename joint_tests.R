rm(list = ls())

library(dplyr)
library(readr)
library(GA)

source("variance_twostage.R")
source("cost_twostage.R")
source("check_constraints.R")
source("repair_joint.R")
source("sample_size_feasibility.R")
source("fitness_ch_kh.R")
source("prepare_strata.R")
source("allocate_twostage_joint.R")

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

total_n <- 900
max_cost <- 20000

ga_result <- allocate_twostage_ga(
  strata = strata,
  total_n = total_n,
  max_cost = max_cost,
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
  allocate_twostage_ga(
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