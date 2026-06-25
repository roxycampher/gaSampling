rm(list = ls())
source("variance_twostage.R")
source("cost_twostage.R")
source("fitness.R")
source("check_constraints.R")
source("prepare_strata.R")

strata <- data.frame(
  stratum = 1:4,
  W_h = c(0.25, 0.25, 0.25, 0.25),
  sigma2_bh = c(10, 8, 12, 6),
  sigma2_wh = c(40, 30, 50, 20),
  cost_cluster = c(100, 120, 90, 110),
  cost_unit = c(10, 8, 12, 9)
)

# example allocation
c_h <- c(10, 12, 8, 5)
k_h <- c(10, 10, 10, 10)

# The stratum sample sizes are
n_h <- c_h * k_h
n_h

(sum(n_h))

# test variance function
var_result <- variance_twostage(
  c_h = c_h,
  k_h = k_h,
  W_h = strata$W_h,
  sigma2_bh = strata$sigma2_bh,
  sigma2_wh = strata$sigma2_wh
)

var_result

# test cost function
cost_result <- cost_twostage(
  c_h = c_h,
  k_h = k_h,
  cost_cluster = strata$cost_cluster,
  cost_unit = strata$cost_unit
)

cost_result

# Test for the constraint function
constraint_result <- check_constraints(
  c_h = c_h,
  k_h = k_h,
  total_n = 350,
  max_cost = 10000,
  cost_cluster = strata$cost_cluster,
  cost_unit = strata$cost_unit
)

constraint_result

# Budget set small to check
check_constraints(
  c_h = c_h,
  k_h = k_h,
  total_n = 350,
  max_cost = 1000,
  cost_cluster = strata$cost_cluster,
  cost_unit = strata$cost_unit
)


# Fixed number of clusters per stratum
c_h <- c(10, 12, 8, 5)

# Candidate k_h values chosen by GA
x_valid <- c(10, 10, 10, 10)

fitness_kh_constrained(
  x = x_valid,
  strata = strata,
  total_n = 350,
  max_cost = 10000,
  c_h = c_h
)

# Testing invalid sample size
x_bad <- c(10, 10, 10, 9)

fitness_kh_constrained(
  x = x_bad,
  strata = strata,
  total_n = 350,
  max_cost = 10000,
  c_h = c_h
)

################################## Data formatting and testing #######################################
library(dplyr)

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

strata

# fitness function

test_alloc <- strata |>
  mutate(
    m_h = 5,   # sample 5 clusters per province
    n_h = 10   # sample 10 units per cluster
  )

fitness_kh_constrained(
  x = x_valid,
  strata = test_alloc,
  total_n = 350,
  max_cost = 10000,
  c_h = c_h
)


# variance


# cost