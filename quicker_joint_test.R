# Here we will use a simulated dataset to test the joint allocation (faster)
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

strata <- data.frame(
  stratum = 1:4,
  C_h = c(25, 20, 30, 18),
  W_h = c(0.25, 0.25, 0.25, 0.25),
  sigma2_bh = c(10, 8, 12, 6),
  sigma2_wh = c(40, 30, 50, 20),
  cost_cluster = c(100, 120, 90, 110),
  cost_unit = c(10, 8, 12, 9)
)

total_n <- 350
max_cost <- 10000

L <- nrow(strata)

initial_c_h <- c(10, 12, 8, 5)
initial_k_h <- c(10, 10, 10, 10)

repaired <- repair_ckh(
  c_h = initial_c_h,
  k_h = initial_k_h,
  total_n = total_n,
  lower_c = rep(1L, L),
  upper_c = strata$C_h,
  lower_k = rep(1L, L),
  upper_k = rep(total_n, L)
)

repaired

fitness_result <- fitness_ckh_constrained(
  x = c(repaired$c_h, repaired$k_h),
  strata = strata,
  total_n = total_n,
  max_cost = max_cost,
  max_c_h = strata$C_h,
  max_k_h = rep(total_n, L)
)

fitness_result

# Test the joint allocation

ga_result <- allocate_twostage_ga_joint(
  strata = strata,
  total_n = total_n,
  max_cost = max_cost,
  popSize = 100,
  maxiter = 200,
  run = 20,
  seed = 123,
  suggestion = c(repaired$c_h, repaired$k_h)
)

ga_result$allocation

ga_result$variance # Check here

ga_result$total_n
ga_result$total_cost
ga_result$feasible
ga_result$constraint_check