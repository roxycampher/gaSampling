rm(list = ls())

library(dplyr)
library(readr)
library(R2BEAT)

source("variance_twostage.R")
source("allocate_twostage_joint.R")
source("fitness_ch_kh.R")
source("cost_twostage.R")
source("read_survey_data.R")
source("sampling_frame.R")
source("repair_joint.R")
source("check_constraints.R")

################################################################################
# 1. PREPARE FRAME & R2BEAT DATA (unchanged — identical to main script)
################################################################################

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

r2b_input <- prepareInputToAllocation1(
  samp_frame  = r2b_frame,
  id_PSU      = "R2B_PSU",
  id_SSU      = "R2B_SSU",
  strata_var  = "R2B_STRATUM",
  target_vars = "count_children_present",
  deff_var    = "R2B_STRATUM",
  domain_var  = "R2B_DOMAIN",
  minimum     = minimum_r2b,
  delta       = 1,
  deff_sugg   = 1
)

run_r2beat <- function(cv_target, input, min_PSU_strat = 2) {
  errors <- data.frame(DOM = "DOM1", CV1 = cv_target)
  beat.2st(
    stratif     = input$strata,
    errors      = errors,
    des_file    = input$des_file,
    psu_file    = input$psu_file,
    rho         = input$rho,
    deft_start  = NULL,
    effst       = input$effst,
    minPSUstrat = min_PSU_strat,
    minnumstrat = 2
  )
}

eval_r2beat_cost <- function(cv_target, input, strata_df, province_lkp, min_r2b) {
  fit <- try(run_r2beat(cv_target = cv_target, input = input), silent = TRUE)
  if (inherits(fit, "try-error")) return(list(cost = NA, fit = NULL, alloc = NULL))

  alloc_tab <- fit$alloc[fit$alloc$STRATUM != "Total", ]
  r2b_nh <- alloc_tab |>
    transmute(R2B_STRATUM = as.integer(as.character(STRATUM)), n_h_r2beat = as.numeric(ALLOC))

  r2b_ch <- lapply(
    split(fit$psu_trs, fit$psu_trs$STRATUM),
    function(x) {
      n_SR  <- sum(x$AR == 1)
      n_NSR <- length(unique(x$SUBSTRAT[x$AR == 0])) * min_r2b
      data.frame(R2B_STRATUM = as.integer(as.character(x$STRATUM[1])), c_h_r2beat = as.numeric(n_SR + n_NSR))
    }
  ) |> bind_rows()

  comb <- r2b_nh |>
    left_join(r2b_ch, by = "R2B_STRATUM") |>
    left_join(province_lkp, by = "R2B_STRATUM") |>
    arrange(match(province, strata_df$province)) |>
    mutate(k_h_r2beat = n_h_r2beat / c_h_r2beat)

  cost <- cost_twostage(
    c_h          = comb$c_h_r2beat,
    k_h          = comb$k_h_r2beat,
    cost_cluster = strata_df$cost_cluster,
    cost_unit    = strata_df$cost_unit
  )

  list(cost = cost, fit = fit, alloc = comb)
}

calibrate_r2beat_bisection <- function(target_cost, input, strata_df, province_lkp, min_r2b,
                                       cv_low = 0.001, cv_high = 0.50, max_iter = 20, tol = 0.01) {
  for (iter in 1:max_iter) {
    cv_mid <- (cv_low + cv_high) / 2
    res_mid <- eval_r2beat_cost(cv_mid, input, strata_df, province_lkp, min_r2b)
    if (is.na(res_mid$cost)) break
    cost_diff <- res_mid$cost - target_cost
    if (abs(cost_diff) / target_cost < tol || (cv_high - cv_low) < 0.0001) {
      cat(sprintf("Bisection converged at Iteration %d! CV: %.5f | Realized Cost: %.0f | Target Cost: %.0f\n",
                  iter, cv_mid, res_mid$cost, target_cost))
      return(res_mid)
    }
    if (res_mid$cost > target_cost) cv_low <- cv_mid else cv_high <- cv_mid
  }
  return(eval_r2beat_cost((cv_low + cv_high) / 2, input, strata_df, province_lkp, min_r2b))
}

integerize_r2beat_allocation <- function(r2beat_alloc, strata_df, target_n,
                                          max_cost, cost_cluster_val, cost_unit_val) {
  L <- nrow(r2beat_alloc)
  lower_c <- rep(1L, L)
  upper_c <- strata_df$C_h
  lower_k <- rep(1L, L)
  upper_k <- rep(target_n, L)

  repair_ckh(
    c_h              = r2beat_alloc$c_h_r2beat,
    k_h              = r2beat_alloc$k_h_r2beat,
    total_n          = target_n,
    lower_c          = lower_c,
    upper_c          = upper_c,
    lower_k          = lower_k,
    upper_k          = upper_k,
    max_iter         = 1000,
    k_tolerance      = 5L,
    c_search_radius  = 2L,
    max_cost         = max_cost,
    cost_cluster     = cost_cluster_val,
    cost_unit        = cost_unit_val
  )
}

################################################################################
# 2. SINGLE-TIER SMOKE TEST
################################################################################

survey_data <- read_survey_data(
  data       = toy_data,
  stratum_id = "province",
  cluster_id = "id_ward",
  y_var      = "count_children_present",
  na_action  = "warn"
)

# --- pick ONE tier from the real script's vectors (index 1, smallest budget) ---
current_max_cost <- 500000
target_n         <- 1000
cost_cluster_val <- 500
cost_unit_val    <- 50

available_c  <- nrow(unique(survey_data[c("strata", "cluster")]))
affordable_c <- floor((current_max_cost - cost_unit_val * target_n) / cost_cluster_val)
full_frame_c <- min(available_c, pmax(1, affordable_c))

full_frame <- prepare_sampling_frame(
  data         = survey_data,
  y_var        = "count_children_present",
  total_c      = full_frame_c,
  cost_cluster = cost_cluster_val,
  cost_unit    = cost_unit_val
)
strata <- full_frame$strata

cat(sprintf("\n=== SMOKE TEST | Budget: $%.0f | target_n: %d ===\n", current_max_cost, target_n))

# --- GA with tiny settings, just to confirm the pipeline runs end to end ---
ga_result <- allocate_twostage_ga_joint(
  strata   = strata,
  total_n  = target_n,
  max_cost = current_max_cost,
  popSize  = 20,     # was 100
  maxiter  = 30,      # was 500
  run      = 5,       # was 30
  seed     = 123
)

cat(sprintf("GA Realized Cost: $%.0f | GA Variance: %.4f\n", ga_result$total_cost, ga_result$variance))

r2beat_matched <- calibrate_r2beat_bisection(
  target_cost  = ga_result$total_cost,
  input        = r2b_input,
  strata_df    = strata,
  province_lkp = province_lookup,
  min_r2b      = minimum_r2b
)

r2beat_integer <- integerize_r2beat_allocation(
  r2beat_alloc     = r2beat_matched$alloc,
  strata_df        = strata,
  target_n         = target_n,
  max_cost         = ga_result$total_cost,
  cost_cluster_val = cost_cluster_val,
  cost_unit_val    = cost_unit_val
)

var_r2beat_continuous <- variance_twostage(
  c_h = r2beat_matched$alloc$c_h_r2beat, k_h = r2beat_matched$alloc$k_h_r2beat,
  W_h = strata$W_h, sigma2_bh = strata$sigma2_bh, sigma2_wh = strata$sigma2_wh
)

var_r2beat_integer <- variance_twostage(
  c_h = r2beat_integer$c_h, k_h = r2beat_integer$k_h,
  W_h = strata$W_h, sigma2_bh = strata$sigma2_bh, sigma2_wh = strata$sigma2_wh
)

test_result <- data.frame(
  GA_Cost = ga_result$total_cost,
  R2BEAT_Cost_Cont = r2beat_matched$cost,
  R2BEAT_Cost_Integer = cost_twostage(r2beat_integer$c_h, r2beat_integer$k_h, cost_cluster_val, cost_unit_val),
  GA_Variance = ga_result$variance,
  R2BEAT_Variance_Cont = var_r2beat_continuous,
  R2BEAT_Variance_Integer = var_r2beat_integer
)

print(test_result)