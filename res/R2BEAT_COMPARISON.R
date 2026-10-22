################################################################################
# GA vs R2BEAT: cost-matched variance comparison across budget tiers
# skips tiers that already have a saved checkpoint.
################################################################################

rm(list = ls())

library(dplyr)
library(readr)
library(R2BEAT)
library(ggplot2)
library(parallel)
library(scales)

source("variance_twostage.R")
source("allocate_twostage_joint.R")
source("fitness_ch_kh.R")
source("cost_twostage.R")
source("read_survey_data.R")
source("sampling_frame.R")
source("repair_joint.R")
source("check_constraints.R")
source("feasible_joint_start.R")


################################################################################
# 0. OUTPUT PATHS
################################################################################

results_dir <- "local_results"
checkpoint_dir <- file.path(results_dir, "checkpoints")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

table_path <- file.path(results_dir, "ga_vs_r2beat_comparison.csv")
plot_variance_path <- file.path(results_dir, "variance_comparison.png")
plot_variance_lines_path <- file.path(results_dir, "variance_by_seed_lines.png")
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
  upper_k <- strata_df$max_k_h   # same realistic per-stratum cap the GA now uses —
  # keeps this a fair comparison rather than letting
  # R2BEAT's integerized allocation propose k_h values
  # no ward can actually supply
  
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

run_one_tier <- function(b, seed, current_max_cost, target_n, cost_cluster_val, cost_unit_val,
                         survey_data, r2b_input, province_lookup, minimum_r2b,
                         ga_popSize = 100, ga_maxiter = 500, ga_run = 30,
                         use_neyman_start = TRUE) {
  
  cat(sprintf("\n======================================================\n"))
  cat(sprintf(" TIER %d | SEED %d | Budget: R%.0f | target_n: %d | C_cluster: R%.0f | C_unit: R%.0f\n",
              b, seed, current_max_cost, target_n, cost_cluster_val, cost_unit_val))
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
  # ---------------------------------------------------------------------------
  # GA STARTING SUGGESTION: NEYMAN c_h ALLOCATION
  # ---------------------------------------------------------------------------
  # The joint GA can run without a user-supplied suggestion, in which case
  # allocate_twostage_ga_joint() builds a simple internal starting point. For
  # this R2BEAT budget-tier comparison, we want the GA to start from a more
  # principled and already-repaired allocation.
  #
  # find_feasible_joint_start() constructs several possible c_h starting shapes,
  # repairs each one into paired c_h/k_h values, and checks feasibility against
  # this tier's resource envelope:
  #
  #   - total sample size close to target_n
  #   - total cost close to current_max_cost
  #   - c_h and k_h inside their stratum bounds
  #
  # The default_starting_suggestion below selects "3. Neyman proportional".
  # ---------------------------------------------------------------------------
  
  no_neyman_start <- function() {
    list(
      neyman_start = list(starting_suggestion = NA_character_, total_c = NA_integer_),
      neyman_start_diag = data.frame(
        repaired_total_n = NA_integer_,
        repaired_total_cost = NA_real_,
        sample_size_gap = NA_integer_,
        feasible = NA
      ),
      ga_suggestion = NULL
    )
  }
  
  if (use_neyman_start) {
    neyman_attempt <- tryCatch(
      {
        ns <- find_feasible_joint_start(
          strata = strata,
          total_n = target_n,
          max_cost = current_max_cost,
          max_candidate_c = full_frame_c,
          default_starting_suggestion = "3. Neyman proportional",
          require_all_starts = FALSE
        )
        
        ns_diag <- ns$diagnostics |>
          filter(starting_suggestion == "3. Neyman proportional")
        
        cat(sprintf(
          "GA Neyman start: c = %d, n = %d, cost = %.0f, sample gap = %d, feasible = %s\n",
          ns$total_c,
          ns_diag$repaired_total_n,
          ns_diag$repaired_total_cost,
          ns_diag$sample_size_gap,
          ns_diag$feasible
        ))
        
        list(neyman_start = ns, neyman_start_diag = ns_diag, ga_suggestion = ns$suggestion)
      },
      error = function(e) {
        # find_feasible_joint_start() can error outright (rather than return
        # an infeasible result) when no repaired shape fits this tier's
        # budget/target_n/cost combination. Rather than losing the whole
        # tier, fall back to letting the GA build its own internal start —
        # same fallback used when use_neyman_start = FALSE.
        warning(sprintf(
          "Tier %d, seed %d: find_feasible_joint_start() failed (%s) — falling back to GA default start.",
          b, seed, conditionMessage(e)
        ))
        no_neyman_start()
      }
    )
    neyman_start <- neyman_attempt$neyman_start
    neyman_start_diag <- neyman_attempt$neyman_start_diag
    ga_suggestion <- neyman_attempt$ga_suggestion
  } else {
    # Quick-test mode: skip find_feasible_joint_start() entirely and let
    # allocate_twostage_ga_joint() build its own internal starting point.
    cat("Skipping Neyman starting suggestion (use_neyman_start = FALSE).\n")
    fallback <- no_neyman_start()
    neyman_start <- fallback$neyman_start
    neyman_start_diag <- fallback$neyman_start_diag
    ga_suggestion <- fallback$ga_suggestion
  }
  
  ga_result <- allocate_twostage_ga_joint(
    strata   = strata,
    total_n  = target_n,
    max_cost = current_max_cost,
    popSize  = ga_popSize,
    maxiter  = ga_maxiter,
    run      = ga_run,
    seed     = seed,
    
    # This is where the Neyman helper output enters the GA (NULL when
    # use_neyman_start = FALSE, letting the GA build its own start).
    suggestion = ga_suggestion
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
  
  cat(sprintf("GA Realized Cost: R%.0f | GA Variance: %.4f\n", ga_result$total_cost, ga_result$variance))
  
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
    Seed                    = seed,
    Budget_Target           = current_max_cost,
    C_Cluster               = cost_cluster_val,
    C_Unit                  = cost_unit_val,
    GA_Cost                 = ga_result$total_cost,
    R2BEAT_Cost_Cont        = r2beat_matched$cost,
    R2BEAT_Cost_Integer     = cost_r2beat_integer,
    GA_Starting_Suggestion  = neyman_start$starting_suggestion,
    GA_Start_PSU            = neyman_start$total_c,
    GA_Start_n              = neyman_start_diag$repaired_total_n,
    GA_Start_Cost           = neyman_start_diag$repaired_total_cost,
    GA_Start_Sample_Gap     = neyman_start_diag$sample_size_gap,
    GA_Start_Feasible       = neyman_start_diag$feasible,
    GA_PSU                  = sum(ga_result$allocation$c_h),
    R2BEAT_PSU_Cont         = sum(r2beat_matched$alloc$c_h_r2beat),
    R2BEAT_PSU_Integer      = sum(r2beat_integer$c_h),
    GA_SSU                  = ga_result$total_n,
    R2BEAT_SSU_Cont         = sum(r2beat_matched$alloc$n_h_r2beat),
    R2BEAT_SSU_Integer      = sum(r2beat_integer$n_h),
    GA_Variance             = ga_result$variance,
    GA_Feasible             = isTRUE(ga_result$feasible),
    R2BEAT_Variance_Cont    = var_r2beat_continuous,
    R2BEAT_Variance_Integer = var_r2beat_integer
  )
}

################################################################################
# 3. TIER DEFINITIONS
################################################################################

target_costs      <- c(1000000, 3000000, 5000000, 7000000, 10000000)
target_n_vec      <- c(   2500,   10000,   17500,   25000,    35000)
cost_cluster_vec  <- c(   1000,    1500,    2000,    2500,     3000)
cost_unit_vec     <- c(    100,     125,     150,     175,      200)

seeds <- c(101, 202, 303, 404, 505)   # 5 seeds per tier

################################################################################
# 4. RUN TIERS WITH CHECKPOINTING
################################################################################

run_and_checkpoint_tier <- function(b, s) {
  checkpoint_file <- file.path(checkpoint_dir, sprintf("tier_%02d_seed_%d.rds", b, s))
  
  if (file.exists(checkpoint_file)) {
    cat(sprintf("Tier %d, seed %d already checkpointed — skipping.\n", b, s))
    return(readRDS(checkpoint_file))
  }
  
  tier_result <- tryCatch(
    run_one_tier(
      b                 = b,
      seed              = s,
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
      warning(sprintf("Tier %d, seed %d failed with error: %s", b, s, conditionMessage(e)))
      NULL
    }
  )
  
  if (!is.null(tier_result)) {
    saveRDS(tier_result, checkpoint_file)
    cat(sprintf("Tier %d, seed %d checkpointed to %s\n", b, s, checkpoint_file))
  }
  
  tier_result
}

tier_seed_grid <- expand.grid(
  b = seq_along(target_costs),
  s = seeds
)

# Number of parallel workers: respect the Slurm allocation if present,
# otherwise fall back to (detected cores - 1), minimum 1.
n_cores <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA)))
if (is.na(n_cores)) {
  n_cores <- max(1L, parallel::detectCores() - 1L)
}
n_cores <- min(n_cores, nrow(tier_seed_grid))
cat(sprintf("Running %d tier x seed combinations across %d core(s).\n", nrow(tier_seed_grid), n_cores))

if (.Platform$OS.type == "unix") {
  all_tier_results <- parallel::mclapply(
    seq_len(nrow(tier_seed_grid)),
    function(i) run_and_checkpoint_tier(tier_seed_grid$b[i], tier_seed_grid$s[i]),
    mc.cores = n_cores,
    mc.preschedule = FALSE
  )
} else {
  # mclapply forks and is Linux/macOS-only; fall back to a serial,
  # still-checkpointed loop if ever run on Windows.
  all_tier_results <- lapply(
    seq_len(nrow(tier_seed_grid)),
    function(i) run_and_checkpoint_tier(tier_seed_grid$b[i], tier_seed_grid$s[i])
  )
}

# mclapply returns a try-error object for a crashed worker instead of
# stopping the whole job; detect and retry those serially. Already-
# checkpointed combinations from other workers are simply reloaded.
failed <- vapply(all_tier_results, function(x) inherits(x, "try-error"), logical(1))
if (any(failed)) {
  cat(sprintf("%d combination(s) failed in parallel and will be retried serially.\n", sum(failed)))
  for (i in which(failed)) {
    all_tier_results[[i]] <- run_and_checkpoint_tier(tier_seed_grid$b[i], tier_seed_grid$s[i])
  }
}

################################################################################
# 5. COMBINE CHECKPOINTS INTO FINAL TABLE
################################################################################

checkpoint_files <- list.files(checkpoint_dir, pattern = "^tier_.*_seed_.*\\.rds$", full.names = TRUE)

if (length(checkpoint_files) == 0) {
  stop("No tier/seed combinations completed successfully — nothing to combine.")
}

final_results <- bind_rows(lapply(checkpoint_files, readRDS)) |>
  arrange(Tier, Seed)

write_csv(final_results, table_path)
cat(sprintf("\nFinal table written to %s\n", table_path))
print(final_results)

################################################################################
# 6. PLOTS
################################################################################

# --- 6a. Variance comparison across budget tiers, spread across seeds ---
variance_long <- final_results |>
  select(Tier, Budget_Target, Seed, GA_Variance, R2BEAT_Variance_Cont, R2BEAT_Variance_Integer) |>
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
    ),
    Budget_Label = factor(
      Budget_Target,
      levels = sort(unique(Budget_Target)),
      labels = scales::label_number(prefix = "R", big.mark = ",")(sort(unique(Budget_Target)))
    )
  )

p_variance <- ggplot(variance_long, aes(x = Budget_Label, y = Variance, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.75), width = 0.6, outlier.size = 1) +
  labs(
    title = "Variance at matched cost: GA vs R2BEAT, across seeds",
    subtitle = sprintf("%d seeds per tier — GA variance only meaningful where GA_Feasible is TRUE", length(unique(final_results$Seed))),
    x = "Budget tier (Rand, matched across methods within tier)",
    y = "Total variance (var_total_h summed)"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

ggsave(plot_variance_path, p_variance, width = 9, height = 5.5, dpi = 300)
cat(sprintf("Variance plot written to %s\n", plot_variance_path))

# --- 6a2. Variance by seed across budget tiers (line plot) ---
# Thin lines trace each individual seed's trajectory across tiers; bold
# lines show the mean per method. Complements the boxplot above: the
# boxplot shows spread at a given budget, this shows whether any single
# seed is consistently high/low across the whole tier range, or whether
# seeds cross over unpredictably (a stability concern the boxplot alone
# wouldn't reveal).
p_variance_lines <- ggplot(
  variance_long,
  aes(x = Budget_Target, y = Variance, color = Method, group = interaction(Method, Seed))
) +
  geom_line(alpha = 0.35, linewidth = 0.6) +
  geom_point(alpha = 0.5, size = 1.5) +
  stat_summary(aes(group = Method), fun = mean, geom = "line", linewidth = 1.4) +
  stat_summary(aes(group = Method), fun = mean, geom = "point", size = 3) +
  scale_x_continuous(labels = scales::label_number(prefix = "R", big.mark = ",")) +
  labs(
    title = "Variance by seed across budget tiers",
    subtitle = "Thin lines = individual seeds; bold lines = mean across seeds per method",
    x = "Budget tier (Rand)",
    y = "Total variance (var_total_h summed)"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

ggsave(plot_variance_lines_path, p_variance_lines, width = 9, height = 5.5, dpi = 300)
cat(sprintf("Per-seed variance line plot written to %s\n", plot_variance_lines_path))

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
  scale_x_continuous(labels = scales::label_number(prefix = "R", big.mark = ",")) +
  scale_y_continuous(labels = scales::label_number(prefix = "R", big.mark = ",")) +
  labs(
    title = "Cost-matching sanity check",
    subtitle = "Points should sit on the dashed line:confirms the bisection actually matched budgets",
    x = "Target budget for tier (Rand)",
    y = "Realized cost (Rand)"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

ggsave(plot_cost_path, p_cost, width = 8, height = 5.5, dpi = 300)
cat(sprintf("Cost-match check plot written to %s\n", plot_cost_path))

cat("\nDone.\n")