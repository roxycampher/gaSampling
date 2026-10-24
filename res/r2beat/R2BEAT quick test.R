rm(list = ls())

library(dplyr)
library(readr)
library(R2BEAT)
library(ggplot2)

source("variance_twostage.R")
source("allocate_twostage_joint.R")
source("fitness_ch_kh.R")
source("cost_twostage.R")
source("read_survey_data.R")
source("sampling_frame.R")
source("repair_joint.R")
source("check_constraints.R")

toy_data <- read.csv("ecdc-2021-public-v2.csv")
toy_data$count_children_present <- readr::parse_number(toy_data$count_children_present)

r2b_frame <- toy_data |>
  filter(!is.na(province), !is.na(id_ward), !is.na(count_children_present))

province_lookup <- r2b_frame |>
  distinct(province) |>
  arrange(province) |>
  mutate(R2B_STRATUM = row_number())

r2b_frame <- r2b_frame |>
  left_join(province_lookup, by = "province") |>
  mutate(
    PSU_KEY    = interaction(province, id_ward, drop = TRUE),
    R2B_PSU    = as.integer(PSU_KEY),
    R2B_SSU    = row_number(),
    R2B_DOMAIN = 1L,
    one        = 1
  )

minimum_r2b <- 1

alloc <- beat.2st(
  stratif = inp$strata,
  errors = cv,
  des_file = inp$des_file,
  psu_file = inp$psu_file,
  rho = inp$rho,
  deft_start = NULL
)

