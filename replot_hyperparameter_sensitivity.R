# Replot saved sensitivity CSVs. This script does not run the GA.
library(dplyr)
library(ggplot2)

# ------------------------------------------------------------
# Folders and appearance
# ------------------------------------------------------------
input_dir <- "C:/Users/roxyc/OneDrive - University of Cape Town/Documents/UCT/2026/Project/gaSampling/hyperparameter_sensitivity_results"
output_dir <- "C:/Users/roxyc/OneDrive - University of Cape Town/Documents/UCT/2026/Project/gaSampling/hyperparameter_faceted_plots"
plot_dpi <- 200

# ------------------------------------------------------------
# Load existing results
# ------------------------------------------------------------
read_results <- function(filename, columns) {
  path <- file.path(input_dir, filename)
  if (!file.exists(path)) stop("Missing CSV: ", path)
  data <- read.csv(path, stringsAsFactors = FALSE)
  missing <- setdiff(columns, names(data))
  if (length(missing)) stop(filename, " is missing: ", paste(missing, collapse = ", "))
  data
}
run_summary <- read_results("run_summary.csv",
                            c("config_id", "phase", "feasible", "variance", "generations", "popSize", "maxiter", "run"))
configuration_summary <- read_results("configuration_summary.csv",
                                      c("config_id", "phase", "feasibility_rate", "fallback_rate", "median_seconds", "median_variance", "maxiter"))
paired <- read_results("paired_comparisons.csv",
                       c("config_id", "phase", "starting_suggestion", "variance_change_pct"))
allocation_results <- read_results("allocations.csv",
                                   c("config_id", "phase", "starting_suggestion", "stratum", "feasible", "c_h", "k_h", "n_h"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Plot variables and shared colours
# ------------------------------------------------------------
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
variance_plot <- boxplot_by_population(filter(run_summary, feasible), "variance",
                                       "Variance by population size.", "Final variance (feasible runs)")
paired_variance_plot <- boxplot_by_population(paired, "variance_change_pct",
                                              "Paired comparison with baseline.", "Variance change from baseline (%)") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey45") +
  labs(subtitle = "Negative values favour the configuration; feasible pairs only")
# Individual run times, not configuration medians, form the runtime boxes.
runtime_data <- run_summary |>
  filter(is.na(error)) |>
  mutate(runtime_hours = elapsed_seconds / 3600)
runtime_plot <- boxplot_by_population(runtime_data, "runtime_hours",
                                      "Runtime by population size.", "Runtime (hours)")
generation_plot <- boxplot_by_population(run_summary, "generations",
                                         "Generations used.", "Actual generations") +
  geom_point(aes(y = maxiter), position = position_dodge(0.8),
             shape = 4, size = 2, na.rm = TRUE) +
  labs(subtitle = "Crosses mark the maximum allowed generations.")
variance_plot
paired_variance_plot
runtime_plot
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
    ggsave(file.path(output_dir, paste0(name, ".", extension)), plot = plot,
           width = width, height = height, dpi = plot_dpi, bg = "white", limitsize = FALSE)
  }
}
plots <- list(variance_boxplots = variance_plot, paired_variance_change = paired_variance_plot,
              runtime_boxplots = runtime_plot, generation_boxplots = generation_plot,
              feasibility_rates = feasibility_plot, fallback_rates = fallback_plot,
              gauteng_ch_stability = gauteng_ch_plot,
              gauteng_kh_stability = gauteng_kh_plot,
              gauteng_nh_stability = gauteng_nh_plot)
for (name in names(plots)) save_plot(plots[[name]], name)
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
          file.path(output_dir, "iteration_colours.csv"), row.names = FALSE)
message("Faceted plots saved to: ", normalizePath(output_dir))

