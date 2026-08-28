#!/usr/bin/env Rscript
################################################################################
# GA vs R2BEAT: cost-matched variance comparison across budget tiers
# Designed to run unattended on HPC. Checkpoints each tier to disk so a
# crash / walltime kill doesn't lose completed work — rerunning the script
# skips tiers that already have a saved checkpoint.
################################################################################

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

################################################################################
# 0. OUTPUT PATHS
################################################################################

results_dir <- "local_results"
checkpoint_dir <- file.path(results_dir, "checkpoints")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

table_path <- file.path(results_dir, "ga_vs_r2beat_comparison.csv")
plot_variance_path <- file.path(results_dir, "variance_comparison.png")
plot_cost_path <- file.path(results_dir, "cost_match_check.png")

################################################################################
# 1. PREPARE FRAME & R2BEAT DATA
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

survey_data <- read_survey_data(
  data       = toy_data,
  stratum_id = "province",
  cluster_id = "id_ward",
  y_var      = "count_children_present",
  na_action  = "warn"
)

################################################################################
# 2. HELPER FUNCTIONS
################################################################################

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
      cat(sprintf("  Bisection converged at iter %d | CV: %.5f | Realized cost: %.0f | Target: %.0f\n",
                  iter, cv_mid, res_mid$cost, target_cost))
      return(res_mid)
    }
    if (res_mid$cost > target_cost) cv_low <- cv_mid else cv_high <- cv_mid
  }
  cat("  Bisection hit max_iter without converging within tolerance — using best available.\n")
  eval_r2beat_cost((cv_low + cv_high) / 2, input, strata_df, province_lkp, min_r2b)
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

run_one_tier <- function(b, current_max_cost, target_n, cost_cluster_val, cost_unit_val,
                          survey_data, r2b_input, province_lookup, minimum_r2b) {

  cat(sprintf("\n======================================================\n"))
  cat(sprintf(" TIER %d | Budget: $%.0f | target_n: %d | C_cluster: $%.0f | C_unit: $%.0f\n",
              b, current_max_cost, target_n, cost_cluster_val, cost_unit_val))
  cat(sprintf("======================================================\n"))

  if ((cost_unit_val * target_n) >= current_max_cost) {
    warning(sprintf("Tier %d: unit costs exceed budget — skipping.", b))
    return(NULL)
  }

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

  ga_result <- allocate_twostage_ga_joint(
    strata   = strata,
    total_n  = target_n,
    max_cost = current_max_cost,
    popSize  = 20,
    maxiter  = 50,
    run      = 20,
    seed     = 123
  )
  
  # if the allocation overshoots
  if (!isTRUE(ga_result$feasible)) {
    warning(sprintf(
      "Tier %d: GA returned infeasible allocation. Cost = %.0f, budget = %.0f, sample gap = %d. Skipping.",
      b,
      ga_result$total_cost,
      current_max_cost,
      ga_result$repair$sample_size_gap
    ))
    return(NULL)
  }
  
  cat(sprintf("GA Realized Cost: $%.0f | GA Variance: %.4f\n", ga_result$total_cost, ga_result$variance))

  r2beat_matched <- calibrate_r2beat_bisection(
    target_cost  = ga_result$total_cost,
    input        = r2b_input,
    strata_df    = strata,
    province_lkp = province_lookup,
    min_r2b      = minimum_r2b
  )

  if (is.null(r2beat_matched$alloc)) {
    warning(sprintf("Tier %d: R2BEAT bisection failed — skipping.", b))
    return(NULL)
  }

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

  cost_r2beat_integer <- cost_twostage(
    c_h = r2beat_integer$c_h, k_h = r2beat_integer$k_h,
    cost_cluster = cost_cluster_val, cost_unit = cost_unit_val
  )

  data.frame(
    Tier                    = b,
    Budget_Target           = current_max_cost,
    C_Cluster               = cost_cluster_val,
    C_Unit                  = cost_unit_val,
    GA_Cost                 = ga_result$total_cost,
    R2BEAT_Cost_Cont        = r2beat_matched$cost,
    R2BEAT_Cost_Integer     = cost_r2beat_integer,
    GA_PSU                  = sum(ga_result$allocation$c_h),
    R2BEAT_PSU_Cont         = sum(r2beat_matched$alloc$c_h_r2beat),
    R2BEAT_PSU_Integer      = sum(r2beat_integer$c_h),
    GA_SSU                  = ga_result$total_n,
    R2BEAT_SSU_Cont         = sum(r2beat_matched$alloc$n_h_r2beat),
    R2BEAT_SSU_Integer      = sum(r2beat_integer$n_h),
    GA_Variance             = ga_result$variance,
    R2BEAT_Variance_Cont    = var_r2beat_continuous,
    R2BEAT_Variance_Integer = var_r2beat_integer
  )
}

################################################################################
# 3. TIER DEFINITIONS
################################################################################

target_costs      <- c(500000, 1000000, 2000000, 5000000, 10000000)
target_n_vec      <- c(  1000,    2500,    5000,   15000,    35000)
cost_cluster_vec  <- c(   500,    1000,    1500,    2000,     3000)
cost_unit_vec     <- c(    50,     100,     120,     150,      200)

################################################################################
# 4. RUN TIERS WITH CHECKPOINTING
################################################################################

for (b in seq_along(target_costs)) {

  checkpoint_file <- file.path(checkpoint_dir, sprintf("tier_%02d.rds", b))

  if (file.exists(checkpoint_file)) {
    cat(sprintf("Tier %d already checkpointed — skipping.\n", b))
    next
  }

  tier_result <- tryCatch(
    run_one_tier(
      b                 = b,
      current_max_cost  = target_costs[b],
      target_n          = target_n_vec[b],
      cost_cluster_val  = cost_cluster_vec[b],
      cost_unit_val     = cost_unit_vec[b],
      survey_data       = survey_data,
      r2b_input         = r2b_input,
      province_lookup   = province_lookup,
      minimum_r2b       = minimum_r2b
    ),
    error = function(e) {
      warning(sprintf("Tier %d failed with error: %s", b, conditionMessage(e)))
      NULL
    }
  )

  if (!is.null(tier_result)) {
    saveRDS(tier_result, checkpoint_file)
    cat(sprintf("Tier %d checkpointed to %s\n", b, checkpoint_file))
  }
}

################################################################################
# 5. COMBINE CHECKPOINTS INTO FINAL TABLE
################################################################################

checkpoint_files <- list.files(checkpoint_dir, pattern = "^tier_.*\\.rds$", full.names = TRUE)

if (length(checkpoint_files) == 0) {
  stop("No tiers completed successfully — nothing to combine.")
}

final_results <- bind_rows(lapply(checkpoint_files, readRDS)) |>
  arrange(Tier)

write_csv(final_results, table_path)
cat(sprintf("\nFinal table written to %s\n", table_path))
print(final_results)

################################################################################
# 6. PLOTS
################################################################################

# --- 6a. Variance comparison across budget tiers ---
variance_long <- final_results |>
  select(Budget_Target, GA_Variance, R2BEAT_Variance_Cont, R2BEAT_Variance_Integer) |>
  tidyr::pivot_longer(
    cols = c(GA_Variance, R2BEAT_Variance_Cont, R2BEAT_Variance_Integer),
    names_to = "Method",
    values_to = "Variance"
  ) |>
  mutate(
    Method = recode(Method,
      GA_Variance                 = "GA (this package)",
      R2BEAT_Variance_Cont        = "R2BEAT (continuous)",
      R2BEAT_Variance_Integer     = "R2BEAT (integer, cost-aware)"
    )
  )

p_variance <- ggplot(variance_long, aes(x = Budget_Target, y = Variance, color = Method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_x_continuous() +
  labs(
    title = "Variance at matched cost: GA vs R2BEAT",
    subtitle = "Lower is better — all methods evaluated at the same realized budget per tier",
    x = "Budget (matched across methods within tier)",
    y = "Total variance (var_total_h summed)"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

ggsave(plot_variance_path, p_variance, width = 8, height = 5.5, dpi = 300)
cat(sprintf("Variance plot written to %s\n", plot_variance_path))

# --- 6b. Cost-match sanity check (confirms bisection actually matched budgets) ---
cost_long <- final_results |>
  select(Budget_Target, GA_Cost, R2BEAT_Cost_Cont, R2BEAT_Cost_Integer) |>
  tidyr::pivot_longer(
    cols = c(GA_Cost, R2BEAT_Cost_Cont, R2BEAT_Cost_Integer),
    names_to = "Method",
    values_to = "RealizedCost"
  ) |>
  mutate(
    Method = recode(Method,
      GA_Cost             = "GA (this package)",
      R2BEAT_Cost_Cont    = "R2BEAT (continuous)",
      R2BEAT_Cost_Integer = "R2BEAT (integer, cost-aware)"
    )
  )

p_cost <- ggplot(cost_long, aes(x = Budget_Target, y = RealizedCost, color = Method)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(size = 3) +
  scale_x_continuous() +
  scale_y_continuous() +
  labs(
    title = "Cost-matching sanity check",
    subtitle = "Points should sit on the dashed line — confirms the bisection actually matched budgets",
    x = "Target budget for tier",
    y = "Realized cost"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

ggsave(plot_cost_path, p_cost, width = 8, height = 5.5, dpi = 300)
cat(sprintf("Cost-match check plot written to %s\n", plot_cost_path))

cat("\nDone.\n")

