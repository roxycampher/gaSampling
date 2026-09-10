rm(list = ls())

library(dplyr)
library(readr)
library(GA)
library(ggplot2)
library(tidyr)

# ------------------------------------------------------------
# Hyperparameter settings
# ------------------------------------------------------------
# Change the settings here before running the sections below.

joint_ga_pop_size <- 200
joint_ga_maxiter <- 500
joint_ga_run <- 50
joint_ga_seed_count <- 20 # Eight configurations x 20 seeds = 160 GA runs
joint_ga_max_candidate_c <- 1000

population_sizes <- c(50, 100, 200, 400)
maximum_iterations <- c(250, 500, 1000)
stopping_patience <- c(25, 50, 100)

include_interactions <- FALSE # Only change one hyperparamter at a time
validation_configuration <- NULL # e.g. "p400_i500_r50", after inspecting results
quick_check <- FALSE # When set to TRUE runs a small check
output_dir <- "hyperparameter_sensitivity_results"