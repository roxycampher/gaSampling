#!
################################################################################
# TEST 2: R2BEAT-first comparison
#
# For each CV target, R2BEAT determines the resource envelope (n, cost)
# it needs to hit that precision. That same envelope (n as floor, cost as
# cap) is then handed to our GA, to see whether a better allocation SHAPE
# achieves lower variance with the exact resources R2BEAT decided were
# necessary — not an independent precision-fixed test, since the GA still
# only ever sees a fixed cost + n floor. 
################################################################################

rm(list = ls())

library(dplyr)
library(readr)
library(R2BEAT)
library(ggplot2)
library(GA)

source("variance_twostage.R")
source("allocate_twostage_joint.R")
source("fitness_ch_kh.R")
source("cost_twostage.R")
source("read_survey_data.R")
source("sampling_frame.R")
source("repair_joint.R")
source("check_constraints.R")

################################################################################
# 0. OUTPUT PATHS (separate from Test 1's — no collision)
################################################################################

results_dir <- "hpc_results_test2"
checkpoint_dir <- file.path(results_dir, "checkpoints")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

table_path <- file.path(results_dir, "ga_vs_r2beat_test2.csv")
plot_variance_path <- file.path(results_dir, "test2_variance_comparison.png")
plot_cost_path <- file.path(results_dir, "test2_cost_check.png")

################################################################################
# 1. PREPARE FRAME & R2BEAT DATA (identical setup to Test 1)
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

# Fixed cost parameters for this test — NOT swept, unlike Test 1's tiers.
# Only the CV target varies; cost_cluster/cost_unit stay constant so any
# variance difference across the CV grid is attributable to precision level
# alone, not a moving cost structure.
cost_cluster_val <- 1000
cost_unit_val    <- 100

available_c <- nrow(unique(survey_data[c("strata", "cluster")]))

full_frame <- prepare_sampling_frame(
  data         = survey_data,
  y_var        = "count_children_present",
  total_c      = available_c,   # generous ceiling — GA's own max_cost constrains it, not this
  cost_cluster = cost_cluster_val,
  cost_unit    = cost_unit_val
)
strata <- full_frame$strata

################################################################################
# 2. HELPER FUNCTIONS (same as Test 1's eval_r2beat_cost / integerize)
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

run_one_cv <- function(idx, cv_target, r2b_input, strata, province_lookup, minimum_r2b,
                        cost_cluster_val, cost_unit_val) {

  cat(sprintf("\n======================================================\n"))
  cat(sprintf(" CV %d | Target CV: %.4f\n", idx, cv_target))
  cat(sprintf("======================================================\n"))

  r2b_res <- eval_r2beat_cost(cv_target, r2b_input, strata, province_lookup, minimum_r2b)
  if (is.na(r2b_res$cost)) {
    warning(sprintf("CV %d: R2BEAT failed at CV = %.4f — skipping.", idx, cv_target))
    return(NULL)
  }

  target_n_i    <- round(sum(r2b_res$alloc$n_h_r2beat))
  target_cost_i <- r2b_res$cost

  cat(sprintf("R2BEAT: n = %d, cost = %.0f\n", target_n_i, target_cost_i))

  ga_res <- allocate_twostage_ga_joint(
    strata   = strata,
    total_n  = target_n_i,
    max_cost = target_cost_i,     # GA gets exactly the resources R2BEAT decided this CV needs
    popSize  = 100,
    maxiter  = 500,
    run      = 30,
    seed     = 123
  )
  
  if (!isTRUE(ga_res$feasible)) {
    warning(sprintf(
      "CV %d: GA returned infeasible allocation. Cost = %.0f, cap = %.0f, sample gap = %d. Skipping.",
      idx,
      ga_res$total_cost,
      target_cost_i,
      ga_res$repair$sample_size_gap
    ))
    return(NULL)
  }
  
  cat(sprintf("GA: n = %d, cost = %.0f, variance = %.4f\n",
              ga_res$total_n, ga_res$total_cost, ga_res$variance))

  r2b_integer <- integerize_r2beat_allocation(
    r2b_res$alloc, strata, target_n_i, target_cost_i, cost_cluster_val, cost_unit_val
  )

  var_r2b_cont <- variance_twostage(
    c_h = r2b_res$alloc$c_h_r2beat, k_h = r2b_res$alloc$k_h_r2beat,
    W_h = strata$W_h, sigma2_bh = strata$sigma2_bh, sigma2_wh = strata$sigma2_wh
  )
  var_r2b_int <- variance_twostage(
    c_h = r2b_integer$c_h, k_h = r2b_integer$k_h,
    W_h = strata$W_h, sigma2_bh = strata$sigma2_bh, sigma2_wh = strata$sigma2_wh
  )
  cost_r2b_int <- cost_twostage(
    c_h = r2b_integer$c_h, k_h = r2b_integer$k_h,
    cost_cluster = cost_cluster_val, cost_unit = cost_unit_val
  )

  data.frame(
    CV_Target               = cv_target,
    R2BEAT_n                = target_n_i,
    R2BEAT_Cost              = target_cost_i,
    R2BEAT_Cost_Integer      = cost_r2b_int,
    GA_n                     = ga_res$total_n,
    GA_Cost                  = ga_res$total_cost,
    GA_Variance              = ga_res$variance,
    R2BEAT_Variance_Cont     = var_r2b_cont,
    R2BEAT_Variance_Integer  = var_r2b_int
  )
}

################################################################################
# 3. CV GRID
################################################################################

cv_grid <- c(0.005, 0.01, 0.02, 0.05, 0.1)

################################################################################
# 4. RUN WITH CHECKPOINTING
################################################################################

for (idx in seq_along(cv_grid)) {

  checkpoint_file <- file.path(checkpoint_dir, sprintf("cv_%02d.rds", idx))

  if (file.exists(checkpoint_file)) {
    cat(sprintf("CV %d already checkpointed — skipping.\n", idx))
    next
  }

  cv_result <- tryCatch(
    run_one_cv(
      idx               = idx,
      cv_target         = cv_grid[idx],
      r2b_input         = r2b_input,
      strata            = strata,
      province_lookup   = province_lookup,
      minimum_r2b       = minimum_r2b,
      cost_cluster_val  = cost_cluster_val,
      cost_unit_val     = cost_unit_val
    ),
    error = function(e) {
      warning(sprintf("CV %d failed with error: %s", idx, conditionMessage(e)))
      NULL
    }
  )

  if (!is.null(cv_result)) {
    saveRDS(cv_result, checkpoint_file)
    cat(sprintf("CV %d checkpointed to %s\n", idx, checkpoint_file))
  }
}

################################################################################
# 5. COMBINE CHECKPOINTS
################################################################################

checkpoint_files <- list.files(checkpoint_dir, pattern = "^cv_.*\\.rds$", full.names = TRUE)

if (length(checkpoint_files) == 0) {
  stop("No CV points completed successfully — nothing to combine.")
}

final_results <- bind_rows(lapply(checkpoint_files, readRDS)) |>
  arrange(CV_Target)

write_csv(final_results, table_path)
cat(sprintf("\nFinal table written to %s\n", table_path))
print(final_results)

################################################################################
# 6. PLOTS
################################################################################

variance_long <- final_results |>
  select(CV_Target, GA_Variance, R2BEAT_Variance_Cont, R2BEAT_Variance_Integer) |>
  tidyr::pivot_longer(
    cols = c(GA_Variance, R2BEAT_Variance_Cont, R2BEAT_Variance_Integer),
    names_to = "Method", values_to = "Variance"
  ) |>
  mutate(Method = recode(Method,
    GA_Variance = "GA (this package)",
    R2BEAT_Variance_Cont = "R2BEAT (continuous)",
    R2BEAT_Variance_Integer = "R2BEAT (integer, cost-aware)"
  ))

p_variance <- ggplot(variance_long, aes(x = CV_Target, y = Variance, color = Method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  labs(
    title = "Test 2: Variance at R2BEAT's own resource envelope",
    subtitle = "For each R2BEAT CV target, both methods get R2BEAT's resulting (n, cost)",
    x = "R2BEAT CV target",
    y = "Total variance (var_total_h summed)"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

ggsave(plot_variance_path, p_variance, width = 8, height = 5.5, dpi = 300)

# Cost-match check: GA_Cost should sit at or just under R2BEAT_Cost by construction
cost_check <- final_results |>
  select(CV_Target, R2BEAT_Cost, GA_Cost)

p_cost <- ggplot(cost_check, aes(x = R2BEAT_Cost, y = GA_Cost)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(size = 3, color = "steelblue") +
  labs(
    title = "Test 2: GA cost vs. its R2BEAT-set budget cap",
    subtitle = "Points at/near the line confirm the GA used the resources it was given",
    x = "R2BEAT realized cost (GA's max_cost)",
    y = "GA realized cost"
  ) +
  theme_minimal(base_size = 13)

ggsave(plot_cost_path, p_cost, width = 8, height = 5.5, dpi = 300)

cat("\nDone.\n")