library(dplyr)
library(readr)
library(GA)

source("variance_twostage.R")
source("cost_twostage.R")
source("check_constraints.R")
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

k_equal <- rep(10, nrow(strata))

k_unequal <- c(8, 5, 18, 14, 11, 7, 6, 4, 17)

variance_twostage(
  c_h = c_h,
  k_h = k_equal,
  W_h = strata$W_h,
  sigma2_bh = strata$sigma2_bh,
  sigma2_wh = strata$sigma2_wh
)

variance_twostage(
  c_h = c_h,
  k_h = k_unequal,
  W_h = strata$W_h,
  sigma2_bh = strata$sigma2_bh,
  sigma2_wh = strata$sigma2_wh
)

cost_twostage(
  c_h = c_h,
  k_h = k_equal,
  cost_cluster = strata$cost_cluster,
  cost_unit = strata$cost_unit
)

cost_twostage(
  c_h = c_h,
  k_h = k_unequal,
  cost_cluster = strata$cost_cluster,
  cost_unit = strata$cost_unit
)
