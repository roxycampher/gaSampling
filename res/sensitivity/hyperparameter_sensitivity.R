library(dplyr)
library(readr)
library(GA)
library(ggplot2)
library(tidyr)
library(parallel)

# ------------------------------------------------------------
# Hyperparameter settings
# ------------------------------------------------------------

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
output_dir <- "res/hyperparameter_sensitivity_results"

source("R/variance_twostage.R")
source("R/cost_twostage.R")
source("R/check_constraints.R")
source("R/repair_joint.R")
source("R/fitness_ch_kh.R")
source("R/feasible_joint_start.R")
source("R/allocate_twostage_joint.R")
source("R/read_survey_data.R")
source("R/sampling_frame.R")
source("R/frame_c.R")

# These lists identify the code and package versions for saved checkpoints.
required_packages <- c("GA", "dplyr", "readr", "ggplot2", "tidyr")

helper_files <- c("variance_twostage.R", "cost_twostage.R", "check_constraints.R",
                  "repair_joint.R", "fitness_ch_kh.R", "feasible_joint_start.R",
                  "allocate_twostage_joint.R", "read_survey_data.R",
                  "sampling_frame.R", "frame_c.R")

helper_files <- file.path("R", helper_files)

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

raw_data <- read.csv("inst/extdata/ecd_subset.csv")
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
                 script = tools::md5sum("res/sensitivity/hyperparameter_sensitivity.R"))

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
if (!quick_check) {
plot_dir <- file.path(experiment_dir, "faceted_plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
plot_dpi <- 200

###############################################################################
configuration_key <- run_summary |>
  distinct(config_id, popSize, maxiter, run)
if (anyDuplicated(configuration_key$config_id)) stop("Inconsistent configuration IDs")
population_levels <- sort(unique(configuration_key$popSize))
patience_levels <- sort(unique(configuration_key$run))
iteration_levels <- sort(unique(configuration_key$maxiter))
iteration_colours <- c("250" = "#4477AA", "500" = "#EE6677", "1000" = "#228833")
if (!all(as.character(iteration_levels) %in% names(iteration_colours))) {
  iteration_colours <- setNames(grDevices::hcl.colors(length(iteration_levels), "Dark 3"),
                                as.character(iteration_levels))
}
# Join settings by ID so all CSVs can use the same factors.
add_settings <- function(data) {
  if (any(!data$config_id %in% configuration_key$config_id)) stop("Unknown configuration ID")
  data |>
    select(-any_of(c("popSize", "maxiter", "run"))) |>
    left_join(configuration_key, by = "config_id") |>
    mutate(population = factor(popSize, levels = population_levels),
           patience = factor(run, levels = patience_levels),
           iteration_limit = factor(maxiter, levels = iteration_levels))
}
run_summary <- add_settings(run_summary)
configuration_summary <- add_settings(configuration_summary)
paired <- add_settings(paired)
allocation_results <- add_settings(allocation_results)

plot_theme <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        strip.background = element_rect(fill = "#EDF2F7"),
        strip.text = element_text(face = "bold"))
colour_scales <- function() list(
  scale_fill_manual(values = iteration_colours, drop = FALSE, name = "Maximum generations"),
  scale_colour_manual(values = iteration_colours, drop = FALSE, name = "Maximum generations"),
  scale_x_discrete(drop = FALSE)
)
# Keep a common vertical scale. Separate phases and starting strategies when present.
population_facets <- function(data) {
  variables <- c(if (length(unique(data$phase)) > 1) "phase",
                 if ("starting_suggestion" %in% names(data) &&
                     length(unique(data$starting_suggestion)) > 1) "starting_suggestion",
                 "population")
  facet_wrap(as.formula(paste("~", paste(variables, collapse = "+"))),
             labeller = label_both, ncol = 2)
}
# Explicit groups ensure limits and boxplots align under dodging.
boxplot_by_population <- function(data, variable, title, y_label) {
  ggplot(data, aes(patience, .data[[variable]], fill = iteration_limit,
                   colour = iteration_limit, group = interaction(patience, iteration_limit))) +
    geom_boxplot(position = position_dodge(0.8), width = 0.65, alpha = 0.4,
                 na.rm = TRUE) +
    colour_scales() + population_facets(data) + plot_theme +
    labs(title = title, x = "Stopping patience", y = y_label,
         caption = "Absent groups represent untested combinations.")
}

# ------------------------------------------------------------
# Variance, paired changes, runtime and generations
# ------------------------------------------------------------
variance_plot <- boxplot_by_population(filter(run_summary, feasible, popSize != 50), "variance",
                                       "Variance by population size.", "Final variance (feasible runs)")
variance_p50_plot <- boxplot_by_population(filter(run_summary, feasible, popSize == 50), "variance",
                                           "Variance: population size 50.", "Final variance (feasible runs)")
# Split the saved comparisons without changing their original baseline.
make_paired_plot <- function(data) {
  boxplot_by_population(data, "variance_change_pct",
                        "Paired comparison with baseline.", "Variance change from baseline (%)") +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey45") +
    labs(subtitle = "Negative values favour the configuration; feasible pairs only")
}
paired_variance_plot <- make_paired_plot(filter(paired, popSize != 50))
paired_variance_p50_plot <- make_paired_plot(filter(paired, popSize == 50))
# One point per feasible run; run is the stopping-patience limit.
runtime_data <- run_summary |>
  filter(is.na(error), feasible, is.finite(variance), is.finite(elapsed_seconds)) |>
  mutate(runtime_hours = elapsed_seconds / 3600)
patience_colours <- c("25" = "#4477AA", "50" = "#EE6677", "100" = "#228833")
if (!all(as.character(patience_levels) %in% names(patience_colours))) {
  patience_colours <- setNames(grDevices::hcl.colors(length(patience_levels), "Dark 3"),
                               as.character(patience_levels))
}
runtime_scatter <- function(data) {
  ggplot(data, aes(
    runtime_hours, variance,
    colour = patience,
    shape = factor(maxiter)
  )) +
    geom_point(size = 2, alpha = 0.55) +
    scale_colour_manual(
      values = patience_colours, drop = FALSE,
      name = "Stopping patience (run)"
    ) +
    scale_shape_manual(
      values = c("250" = 16, "500" = 17, "1000" = 15),
      name = "Maximum generations"
    ) +
    population_facets(data) + plot_theme +
    labs(
      title = "Runtime versus final variance by population size.",
      subtitle = "Each point represents one feasible GA run.",
      x = "Runtime (hours)",
      y = "Final variance",
      caption = paste(
        "Shapes indicate maximum generation limits;",
        "overlapping runs may coincide."
      )
    )
}
runtime_plot <- runtime_scatter(filter(runtime_data, popSize != 50))
runtime_p50_plot <- runtime_scatter(filter(runtime_data, popSize == 50))

make_generation_plot <- function(data) {
  boxplot_by_population(data, "generations", "Generations used.", "Actual generations") +
    geom_point(aes(y = maxiter), position = position_dodge(0.8),
               shape = 4, size = 2, na.rm = TRUE) +
    labs(subtitle = "Crosses mark the maximum allowed generations.")
}
generation_plot <- make_generation_plot(filter(run_summary, popSize != 50))
generation_p50_plot <- make_generation_plot(filter(run_summary, popSize == 50))
generation_p50_plot
variance_plot
paired_variance_plot
runtime_plot
variance_p50_plot
paired_variance_p50_plot
runtime_p50_plot
generation_plot

# ------------------------------------------------------------
# Feasibility and fallback remain rate plots
# ------------------------------------------------------------
rate_plot <- function(column, title) {
  data <- configuration_summary |>
    mutate(rate = .data[[column]],
           rate_label = ifelse(is.na(rate), "NA", sprintf("%.0f%%", 100 * rate)),
           label_position = ifelse(is.na(rate), 0, rate))
  # Retain a slot for missing rates without treating them as observed zeros.
  ggplot(data, aes(patience, label_position, fill = iteration_limit,
                   colour = iteration_limit, group = interaction(patience, iteration_limit))) +
    geom_col(position = position_dodge(0.8), width = 0.65, alpha = 0.7) +
    geom_point(aes(alpha = !is.na(rate)), position = position_dodge(0.8), size = 2) +
    scale_alpha_manual(values = c("FALSE" = 0, "TRUE" = 1), guide = "none") +
    geom_text(aes(label = rate_label), position = position_dodge(0.8),
              vjust = -0.5, size = 3, colour = "grey20") +
    scale_y_continuous(limits = c(0, 1.12), breaks = seq(0, 1, 0.25),
                       labels = function(x) paste0(100 * x, "%")) +
    colour_scales() + population_facets(data) + plot_theme +
    labs(title = title, x = "Stopping patience", y = "Proportion of runs",
         caption = "NA means unrecorded, not zero. Absent groups were not tested.")
}
feasibility_plot <- rate_plot("feasibility_rate", "Feasibility by population size")
fallback_plot <- rate_plot("fallback_rate", "Fallback by population size")
feasibility_plot
fallback_plot

# ------------------------------------------------------------
# Allocation figures: one population per figure, faceted by stratum
# ------------------------------------------------------------
allocation_plot_data <- allocation_results |>
  mutate(across(c(c_h, k_h, n_h), ~ifelse(feasible, .x, NA_real_)))
allocation_plot <- function(data, variable, y_label, population_size) {
  variables <- c(if (length(unique(data$phase)) > 1) "phase",
                 if (length(unique(data$starting_suggestion)) > 1) "starting_suggestion",
                 "stratum")
  ggplot(data, aes(patience, .data[[variable]], fill = iteration_limit,
                   colour = iteration_limit, group = interaction(patience, iteration_limit))) +
    geom_boxplot(position = position_dodge(0.8), width = 0.65, alpha = 0.4, na.rm = TRUE) +
    facet_wrap(as.formula(paste("~", paste(variables, collapse = "+"))),
               scales = "free_y", ncol = 3) +
    colour_scales() + plot_theme +
    labs(title = paste(y_label, "- population", population_size),
         x = "Stopping patience", y = y_label,
         subtitle = "Feasible runs only; each stratum panel has its own vertical scale.")
}

# ------------------------------------------------------------
# Gauteng: compare population sizes on the same vertical scale
# ------------------------------------------------------------
gauteng_plot_data <- filter(allocation_plot_data, stratum == "Gauteng")
if (nrow(gauteng_plot_data) == 0) stop("No Gauteng allocations found in allocations.csv")

gauteng_ch_plot <- boxplot_by_population(
  gauteng_plot_data, "c_h", "Gauteng: Cluster allocation stability.",
  "Sampled clusters (c_h)") +
  labs(subtitle = "Feasible runs only; population panels share the same vertical scale")
gauteng_kh_plot <- boxplot_by_population(
  gauteng_plot_data, "k_h", "Gauteng: Within-cluster allocation stability.",
  "Units per sampled cluster (k_h)") +
  labs(subtitle = "Feasible runs only; population panels share the same vertical scale")
gauteng_nh_plot <- boxplot_by_population(
  gauteng_plot_data, "n_h", "Gauteng: Sample-size stability.",
  "Sampled units (n_h = c_h * k_h)") +
  labs(subtitle = "Feasible runs only; population panels share the same vertical scale")

gauteng_ch_plot
gauteng_kh_plot
gauteng_nh_plot

# ------------------------------------------------------------
# Save PNG and PDF figures, without rerunning any GA experiments
# ------------------------------------------------------------
save_plot <- function(plot, name, width = 12, height = 8) {
  for (extension in c("png", "pdf")) {
    ggsave(file.path(plot_dir, paste0(name, ".", extension)), plot = plot,
           width = width, height = height, dpi = plot_dpi, bg = "white", limitsize = FALSE)
  }
}
plots <- list(variance_boxplots = variance_plot, paired_variance_change = paired_variance_plot,
              runtime_vs_variance = runtime_plot, generation_boxplots = generation_plot,
              feasibility_rates = feasibility_plot, fallback_rates = fallback_plot,
              gauteng_ch_stability = gauteng_ch_plot,
              gauteng_kh_stability = gauteng_kh_plot,
              gauteng_nh_stability = gauteng_nh_plot)
for (name in names(plots)) save_plot(plots[[name]], name)
save_plot(variance_p50_plot, "variance_boxplots_p50", width = 7, height = 5)
save_plot(paired_variance_p50_plot, "paired_variance_change_p50", width = 7, height = 5)
save_plot(runtime_p50_plot, "runtime_vs_variance_p50", width = 7, height = 5)
save_plot(generation_p50_plot, "generation_boxplots_p50", width = 7, height = 5)
allocation_plots <- list()
for (population_size in population_levels) {
  population_data <- filter(allocation_plot_data, popSize == population_size)
  panel_count <- nrow(distinct(population_data, phase, starting_suggestion, stratum))
  ch_stability_plot <- allocation_plot(population_data, "c_h", "Sampled clusters (c_h)", population_size)
  kh_stability_plot <- allocation_plot(population_data, "k_h", "Units per cluster (k_h)", population_size)
  nh_stability_plot <- allocation_plot(population_data, "n_h", "Sampled units (n_h)", population_size)
  population_plots <- list(ch_stability = ch_stability_plot,
                           kh_stability = kh_stability_plot, nh_stability = nh_stability_plot)
  for (name in names(population_plots)) {
    filename <- paste0(name, "_p", population_size)
    allocation_plots[[filename]] <- population_plots[[name]]
    save_plot(population_plots[[name]], filename, width = 14,
              height = max(9, 3 * ceiling(panel_count / 3)))
  }
}
write.csv(data.frame(maxiter = names(iteration_colours), colour = unname(iteration_colours)),
          file.path(plot_dir, "iteration_colours.csv"), row.names = FALSE)

message("Faceted plots saved to: ", normalizePath(plot_dir))

} else {
  message("Skipping faceted plots for the quick check.")
}

##############################################################################################


capture.output(sessionInfo(), file = file.path(experiment_dir, "session_info.txt"))
print(configuration_summary)
message("Results saved to: ", normalizePath(experiment_dir))