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

