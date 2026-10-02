rm(list = ls())

library(dplyr)
library(readr)
library(GA)
library(ggplot2)
library(tidyr)
library(parallel)

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

source("variance_twostage.R")
source("cost_twostage.R")
source("check_constraints.R")
source("repair_joint.R")
source("fitness_ch_kh.R")
source("feasible_joint_start.R")
source("allocate_twostage_joint.R")
source("read_survey_data.R")
source("sampling_frame.R")
source("frame_c.R")

# These lists identify the code and package versions for saved checkpoints.
required_packages <- c("GA", "dplyr", "readr", "ggplot2", "tidyr")

helper_files <- c("variance_twostage.R", "cost_twostage.R", "check_constraints.R",
                  "repair_joint.R", "fitness_ch_kh.R", "feasible_joint_start.R",
                  "allocate_twostage_joint.R", "read_survey_data.R",
                  "sampling_frame.R", "frame_c.R")

stopifnot(length(joint_ga_seed_count) == 1L, is.finite(joint_ga_seed_count),
          joint_ga_seed_count >= 1,
          joint_ga_seed_count == as.integer(joint_ga_seed_count))

if (quick_check) output_dir <- file.path(output_dir, "quick_check")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Full ECDC data and fixed sampling problem
# ------------------------------------------------------------
total_n <- 5000L
max_cost <- 1400000
unit_cost <- 100
cluster_cost <- 2500

raw_data <- read.csv("ecdc-2021-public-v2.csv")
raw_data$count_children_present <- readr::parse_number(raw_data$count_children_present)
survey_data <- read_survey_data(raw_data, stratum_id = "province",
                                cluster_id = "id_ward", y_var = "count_children_present",
                                na_action = "warn")

frame_c_h <- frame_c(survey_data, max_cost, unit_cost, cluster_cost, total_n)
full_frame <- prepare_sampling_frame(survey_data, y_var = "count_children_present",
                                     total_c = frame_c_h, cost_cluster = cluster_cost,
                                     cost_unit = unit_cost)
strata <- full_frame$strata
strata # Check the strata before running the GA

# ------------------------------------------------------------
# Feasible starting allocations
# ------------------------------------------------------------
L <- nrow(strata)
lower_c <- rep(2L, L)
upper_c <- strata$C_h
lower_k <- rep(1L, L)
upper_k <- rep(total_n, L)

joint_start <- find_feasible_joint_start(
  strata = strata,
  total_n = total_n,
  max_cost = max_cost,
  lower_c = lower_c,
  upper_c = upper_c,
  lower_k = lower_k,
  upper_k = upper_k,
  max_candidate_c = min(joint_ga_max_candidate_c, frame_c_h)
)
joint_starting_suggestions <- joint_start$all_starts
starting_suggestion_diagnostics <- joint_start$diagnostics
joint_starting_suggestions
starting_suggestion_diagnostics

starts <- split(joint_starting_suggestions, joint_starting_suggestions$starting_suggestion)
default_start <- "3. Neyman proportional"

################################################################################
#                    HYPERPARAMETER CONFIGURATIONS
################################################################################
# The first row is the baseline. Each other row changes one setting.
baseline_configuration <- data.frame(
  popSize = joint_ga_pop_size,
  maxiter = joint_ga_maxiter,
  run = joint_ga_run
)

configurations <- bind_rows(
  baseline_configuration,
  data.frame(popSize = population_sizes, maxiter = joint_ga_maxiter, run = joint_ga_run),
  data.frame(popSize = joint_ga_pop_size, maxiter = maximum_iterations, run = joint_ga_run),
  data.frame(popSize = joint_ga_pop_size, maxiter = joint_ga_maxiter, run = stopping_patience)
) |>
  distinct()

# Investigate whether population size and stopping patience interact.(When set to TRUE)
if (include_interactions) {
  configurations <- unique(rbind(configurations, expand.grid(
    popSize = c(100, 200, 400), maxiter = 1000, run = c(25, 50, 100))))
}

if (quick_check) { # Quick run set to TRUE
  configurations <- data.frame(popSize = 10, maxiter = 2, run = 2)
  joint_ga_seed_count <- 1L
}

configurations$config_id <- with(configurations, paste0("p", popSize, "_i", maxiter, "_r", run))
baseline_id <- configurations$config_id[1]
configurations

# Use the same seeds and the same starting suggestion for every configuration.
seeds <- seq_len(joint_ga_seed_count)
grid <- merge(configurations, data.frame(seed = seeds), by = NULL)

grid$starting_suggestion <- default_start
grid$phase <- "screening"
chosen <- validation_configuration

# Optional: compare the chosen setting with baseline on fresh seeds and all starts.
if (!is.null(chosen) && !quick_check) {
  stopifnot(length(chosen) == 1L, chosen %in% configurations$config_id)
  validation <- merge(configurations[configurations$config_id %in% c(baseline_id, chosen), ],
                      expand.grid(seed = 10000L + seeds,
                                  starting_suggestion = names(starts),
                                  stringsAsFactors = FALSE), by = NULL)
  validation$phase <- "validation"
  grid <- bind_rows(grid, validation)
}

grid
nrow(grid)

# ------------------------------------------------------------
# Checkpoint setup
# ------------------------------------------------------------
manifest <- list(strata = strata, total_n = total_n, max_cost = max_cost,
                 starts = starts, R = R.version.string,
                 packages = lapply(required_packages, function(p) as.character(packageVersion(p))),
                 helpers = tools::md5sum(helper_files),
                 script = tools::md5sum("hyperparameter_sensitivity.R"))

manifest_tmp <- tempfile(tmpdir = output_dir, fileext = ".rds")
saveRDS(manifest, manifest_tmp, version = 2)
experiment_id <- unname(tools::md5sum(manifest_tmp))
unlink(manifest_tmp)
experiment_dir <- file.path(output_dir, experiment_id)
checkpoint_dir <- file.path(experiment_dir, "checkpoints")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(manifest, file.path(experiment_dir, "manifest.rds"))
write.csv(grid, file.path(experiment_dir, "run_grid.csv"), row.names = FALSE)

#################################################################################
#                    PARALELL RUNS
################################################################################

number_of_workers <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
if (is.na(number_of_workers)) {
  number_of_workers <- max(1L, min(4L, parallel::detectCores() - 1L))
}
message("Using ", number_of_workers, " workers (mclapply).")

# mclapply's forked workers swallow warnings/cat()/message() output when run
# under R CMD BATCH, so each run appends its own status line to a log file
# rather than relying on console output to know which combinations completed.
progress_log <- file.path(experiment_dir, "progress.log")

################################################################################
#                    MULTIPLE RUN TESTS
################################################################################

ga_runs <- parallel::mclapply(
  seq_len(nrow(grid)),
  function(i) {
    
    tryCatch({
      
      settings <- grid[i, ]
      checkpoint <- file.path(checkpoint_dir, paste0(settings$phase, "_", settings$config_id,
                                                     "_start", match(settings$starting_suggestion, names(starts)),
                                                     "_seed", settings$seed, ".rds"))
      
      if (file.exists(checkpoint)) return(readRDS(checkpoint))
      
      cat(sprintf("[%s] start run %d/%d: %s, seed %s, %s\n",
                  format(Sys.time(), "%H:%M:%S"), i, nrow(grid), settings$config_id,
                  settings$seed, settings$starting_suggestion),
          file = progress_log, append = TRUE)
      
      start <- starts[[settings$starting_suggestion]]
      start_variance <- variance_twostage(start$repaired_c_h, start$repaired_k_h,
                                          strata$W_h, strata$sigma2_bh, strata$sigma2_wh)
      
      begin <- proc.time()[["elapsed"]]
      
      result <- tryCatch(allocate_twostage_ga_joint(
        strata, total_n, max_cost, popSize = settings$popSize,
        maxiter = settings$maxiter, run = settings$run, seed = settings$seed,
        suggestion = c(start$repaired_c_h, start$repaired_k_h), monitor = FALSE),
        error = function(e) e)
      
      elapsed <- proc.time()[["elapsed"]] - begin
      
      if (inherits(result, "error")) {
        return(list(summary = cbind(settings, feasible = FALSE, elapsed_seconds = elapsed,
                                    error = conditionMessage(result)), allocation = NULL))}
      
      # Some result fields (e.g. repair$iterations, fallback_used) can come back
      # length-0 rather than a scalar on trivial or edge-case runs, which makes
      # cbind()/data.frame() error on a row-count mismatch. Coerce to NA instead.
      scalar_or_na <- function(x) if (length(x) == 0) NA else x[1]
      variance_val <- scalar_or_na(result$variance)
      total_n_val <- scalar_or_na(result$total_n)
      
      summary <- cbind(settings, feasible = scalar_or_na(result$feasible),
                       variance = variance_val,
                       start_variance = start_variance,
                       improvement_pct = if (is.na(variance_val)) NA_real_
                       else 100 * (start_variance - variance_val) / start_variance,
                       elapsed_seconds = elapsed, generations = scalar_or_na(result$ga_object@iter),
                       fallback_used = scalar_or_na(result$fallback_used),
                       total_n = total_n_val, total_cost = scalar_or_na(result$total_cost),
                       sample_size_gap = if (is.na(total_n_val)) NA_real_ else total_n_val - total_n,
                       repair_iterations = scalar_or_na(result$repair$iterations), error = NA_character_)
      
      allocation <- cbind(settings[rep(1, nrow(result$allocation)), ], result$allocation,
                          feasible = scalar_or_na(result$feasible))
      
      saved <- list(summary = summary, allocation = allocation)
      
      # Interrupted writes leave only a temporary file, which is never reloaded.
      temporary <- paste0(checkpoint, ".tmp")
      saveRDS(saved, temporary)
      if (!file.rename(temporary, checkpoint)) stop("Cannot save checkpoint: ", checkpoint)
      
      cat(sprintf("[%s] done run %d/%d: %s, seed %s, %s, feasible=%s, elapsed=%.1fs\n",
                  format(Sys.time(), "%H:%M:%S"), i, nrow(grid), settings$config_id,
                  settings$seed, settings$starting_suggestion,
                  isTRUE(saved$summary$feasible[1]), elapsed),
          file = progress_log, append = TRUE)
      saved
      
    }, error = function(e) {
      # Catches failures anywhere in the run (not just inside the GA call
      # itself), so a genuine crash never surfaces as an opaque mclapply
      # try-error with no message.
      cat(sprintf("[%s] CRASHED run %d/%d: %s\n",
                  format(Sys.time(), "%H:%M:%S"), i, nrow(grid), conditionMessage(e)),
          file = progress_log, append = TRUE)
      list(summary = cbind(grid[i, ], feasible = FALSE, elapsed_seconds = NA_real_,
                           error = conditionMessage(e)), allocation = NULL)
    })
  },                          # Ends function(i)
  mc.cores = number_of_workers,
  mc.preschedule = FALSE,
  mc.set.seed = TRUE
)

# Flag any run that mclapply itself failed on (e.g. a worker crash/OOM kill),
# which comes back as a try-error rather than the expected list structure —
# distinct from the tryCatch()-handled GA errors already logged per run above.
mclapply_failures <- which(vapply(ga_runs, function(x) inherits(x, "try-error"), logical(1)))
if (length(mclapply_failures) > 0) {
  warning("mclapply worker failure on run(s): ", paste(mclapply_failures, collapse = ", "),
          " -- rerun the script to retry from checkpoint.")
}

# ------------------------------------------------------------
# Run-level results
# ------------------------------------------------------------
if (length(mclapply_failures) > 0) ga_runs <- ga_runs[-mclapply_failures]

run_summary <- bind_rows(lapply(ga_runs, function(x) x$summary))
run_summary

write.csv(run_summary, file.path(experiment_dir, "run_summary.csv"), row.names = FALSE)
if (!"variance" %in% names(run_summary)) stop("All runs failed; inspect run_summary.csv")
allocation_results <- bind_rows(lapply(ga_runs, function(x) x$allocation))
write.csv(allocation_results, file.path(experiment_dir, "allocations.csv"), row.names = FALSE)

################################################################################
#                    SUMMARY TABLES
################################################################################
# Keep infeasible and failed runs in the feasibility rate. Summarise variance
# only for feasible runs. Missing summaries remain NA if none are feasible.
safe_median <- function(x) if (any(is.finite(x))) median(x[is.finite(x)]) else NA_real_

safe_iqr <- function(x) if (any(is.finite(x))) IQR(x[is.finite(x)]) else NA_real_

configuration_summary <- run_summary |>
  group_by(phase, config_id, popSize, maxiter, run) |>
  summarise(runs = n(), errors = sum(!is.na(error)),
            feasibility_rate = mean(feasible),
            median_variance = safe_median(variance[feasible]),
            iqr_variance = safe_iqr(variance[feasible]),
            median_improvement_pct = safe_median(improvement_pct[feasible]),
            median_seconds = safe_median(elapsed_seconds[is.na(error)]),
            median_generations = safe_median(generations),
            fallback_rate = if (all(is.na(fallback_used))) NA_real_ else mean(fallback_used, na.rm = TRUE),
            .groups = "drop")
write.csv(configuration_summary, file.path(experiment_dir, "configuration_summary.csv"), row.names = FALSE)

configuration_summary

# Allocation stability within each stratum and starting strategy.
allocation_summary <- allocation_results |>
  filter(feasible) |>
  group_by(phase, config_id, starting_suggestion, stratum) |>
  summarise(across(c(c_h, k_h, n_h), list(median = median, iqr = IQR)), .groups = "drop")
write.csv(allocation_summary, file.path(experiment_dir, "allocation_summary.csv"), row.names = FALSE)

allocation_summary

# Same-seed, same-start baseline comparison (only pairs with feasible outputs)
baseline <- run_summary |>
  filter(config_id == baseline_id) |>
  select(phase, seed, starting_suggestion, baseline_variance = variance,
         baseline_feasible = feasible)

paired <- run_summary |>
  left_join(baseline, by = c("phase", "seed", "starting_suggestion")) |>
  mutate(variance_change_pct = ifelse(feasible & baseline_feasible,
                                      100 * (variance / baseline_variance - 1), NA_real_))
write.csv(paired, file.path(experiment_dir, "paired_comparisons.csv"), row.names = FALSE)

paired

################################################################################
#                    DIAGNOSTIC PLOTS
################################################################################

variance_plot <- ggplot2::ggplot(
  filter(run_summary, feasible), ggplot2::aes(config_id, variance)) +
  ggplot2::geom_boxplot() + ggplot2::facet_wrap(~phase, scales = "free_x") +
  ggplot2::coord_flip() + ggplot2::theme_bw() +
  ggplot2::labs(x = "Configuration", y = "Variance (feasible runs only)")

variance_plot

feasibility_plot <- ggplot2::ggplot(configuration_summary,
                                    ggplot2::aes(config_id, feasibility_rate)) + ggplot2::geom_col() +
  ggplot2::facet_wrap(~phase, scales = "free_x") + ggplot2::coord_flip() +
  ggplot2::ylim(0, 1) + ggplot2::theme_bw()

feasibility_plot

runtime_plot <- ggplot2::ggplot(configuration_summary,
                                ggplot2::aes(median_seconds, median_variance, label = config_id)) +
  ggplot2::geom_point() + ggplot2::geom_text(vjust = -0.7, size = 3, check_overlap = TRUE) +
  ggplot2::facet_wrap(~phase) + ggplot2::theme_bw() +
  ggplot2::labs(x = "Median runtime (seconds)", y = "Median feasible variance")

runtime_plot

# Paired comparisons use the same seed and start as the baseline.
# Negative values mean lower variance than baseline, zero means no change.
paired_variance_plot <- ggplot2::ggplot(
  paired, ggplot2::aes(config_id, variance_change_pct)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  ggplot2::geom_boxplot(na.rm = TRUE) +
  ggplot2::geom_point(size = 1, alpha = 0.5, na.rm = TRUE) +
  ggplot2::facet_wrap(~starting_suggestion, scales = "free_x") +
  ggplot2::coord_flip() + ggplot2::theme_bw() +
  ggplot2::labs(x = "Configuration", y = "Variance change from baseline (%)")

paired_variance_plot

# Fallback means the GA returned the starting allocation after its own
# repaired solution failed the feasibility check. Failed runs are excluded.
fallback_plot <- ggplot2::ggplot(configuration_summary,
                                 ggplot2::aes(config_id, fallback_rate)) +
  ggplot2::geom_col(na.rm = TRUE) +
  ggplot2::facet_wrap(~phase, scales = "free_x") +
  ggplot2::coord_flip() + ggplot2::ylim(0, 1) + ggplot2::theme_bw() +
  ggplot2::labs(x = "Configuration", y = "Fallback proportion")

fallback_plot

# Include completed infeasible runs too (they still used computation)
# Red crosses mark the maximum generations allowed for each configuration.
generation_plot <- ggplot2::ggplot(run_summary,
                                   ggplot2::aes(config_id, generations)) +
  ggplot2::geom_boxplot(na.rm = TRUE) +
  ggplot2::geom_point(size = 1, alpha = 0.5, na.rm = TRUE) +
  ggplot2::geom_point(data = configuration_summary,
                      ggplot2::aes(y = maxiter), colour = "firebrick", shape = 4, size = 3) +
  ggplot2::facet_wrap(~phase, scales = "free_x") +
  ggplot2::coord_flip() + ggplot2::theme_bw() +
  ggplot2::labs(x = "Configuration", y = "Actual generations")

generation_plot

# Allocation stability across seeds, separately for each stratum and start.
# Keep all rows for facet labels, but omit infeasible values from the boxes.
allocation_plot_data <- allocation_results |>
  mutate(across(c(c_h, k_h, n_h), ~ifelse(feasible, .x, NA_real_)))

ch_stability_plot <- ggplot2::ggplot(allocation_plot_data,
                                     ggplot2::aes(config_id, c_h)) +
  ggplot2::geom_boxplot(na.rm = TRUE) +
  ggplot2::geom_point(size = 1, alpha = 0.4, na.rm = TRUE) +
  ggplot2::facet_wrap(~starting_suggestion + stratum, scales = "free_y", ncol = 3) +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1)) +
  ggplot2::labs(x = "Configuration", y = "Sampled clusters (c_h)")

ch_stability_plot

kh_stability_plot <- ggplot2::ggplot(allocation_plot_data,
                                     ggplot2::aes(config_id, k_h)) +
  ggplot2::geom_boxplot(na.rm = TRUE) +
  ggplot2::geom_point(size = 1, alpha = 0.4, na.rm = TRUE) +
  ggplot2::facet_wrap(~starting_suggestion + stratum, scales = "free_y", ncol = 3) +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1)) +
  ggplot2::labs(x = "Configuration", y = "Units per sampled cluster (k_h)")

kh_stability_plot

nh_stability_plot <- ggplot2::ggplot(allocation_plot_data,
                                     ggplot2::aes(config_id, n_h)) +
  ggplot2::geom_boxplot(na.rm = TRUE) +
  ggplot2::geom_point(size = 1, alpha = 0.4, na.rm = TRUE) +
  ggplot2::facet_wrap(~starting_suggestion + stratum, scales = "free_y", ncol = 3) +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1)) +
  ggplot2::labs(x = "Configuration", y = "Sampled units (n_h = c_h * k_h)")

nh_stability_plot

# ------------------------------------------------------------
# Save plots and session details
# ------------------------------------------------------------
ggsave(file.path(experiment_dir, "variance_boxplots.png"),
       plot = variance_plot, width = 10, height = 6, dpi = 160)
ggsave(file.path(experiment_dir, "feasibility_rates.png"),
       plot = feasibility_plot, width = 10, height = 6, dpi = 160)
ggsave(file.path(experiment_dir, "runtime_vs_variance.png"),
       plot = runtime_plot, width = 10, height = 6, dpi = 160)
ggsave(file.path(experiment_dir, "paired_variance_change.png"),
       plot = paired_variance_plot, width = 12, height = 8, dpi = 160)
ggsave(file.path(experiment_dir, "fallback_rates.png"),
       plot = fallback_plot, width = 10, height = 6, dpi = 160)
ggsave(file.path(experiment_dir, "generation_boxplots.png"),
       plot = generation_plot, width = 10, height = 6, dpi = 160)
allocation_panel_count <- allocation_plot_data |>
  distinct(phase, starting_suggestion, stratum) |>
  nrow()
allocation_plot_height <- max(9, 3 * ceiling(allocation_panel_count / 3))
ggsave(file.path(experiment_dir, "ch_stability.png"),
       plot = ch_stability_plot, width = 14, height = allocation_plot_height, dpi = 160)
ggsave(file.path(experiment_dir, "kh_stability.png"),
       plot = kh_stability_plot, width = 14, height = allocation_plot_height, dpi = 160)
ggsave(file.path(experiment_dir, "nh_stability.png"),
       plot = nh_stability_plot, width = 14, height = allocation_plot_height, dpi = 160)
capture.output(sessionInfo(), file = file.path(experiment_dir, "session_info.txt"))
print(configuration_summary)
message("Results saved to: ", normalizePath(experiment_dir))