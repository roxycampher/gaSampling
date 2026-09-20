# Replot saved sensitivity CSVs. This script does not run the GA.
library(dplyr)
library(ggplot2)

# ------------------------------------------------------------
# Folders and appearance
# ------------------------------------------------------------
input_dir <- "C:/Users/roxyc/OneDrive - University of Cape Town/Documents/UCT/2026/Project/gaSampling/hyperparameter_sensitivity_results"
output_dir <- "C:/Users/roxyc/OneDrive - University of Cape Town/Documents/UCT/2026/Project/gaSampling/hyperparameter_colour_plots"
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
# One named colour mapping shared by every figure
# ------------------------------------------------------------
# Order numerically by population, iteration limit and stopping patience.
configuration_key <- run_summary |>
  distinct(config_id, popSize, maxiter, run) |>
  arrange(popSize, maxiter, run)
config_levels <- configuration_key$config_id
if (anyDuplicated(config_levels)) stop("A configuration ID has inconsistent settings.")

# Sixteen distinct colours for the current experiment. Edit hex codes here
# to customise them. Larger experiments use an automatically expanded palette.

palette <- c("#4477AA", "#EE6677", "#228833", "#AA3377",
             "#CCBB44", "#66CCEE", "#332288", "#CC6677",
             "#117733", "#882255", "#44AA99", "#DDCC77",
             "#661100", "#6699CC", "#AA4499", "#E17C05")
if (length(config_levels) > length(palette)) {
  palette <- grDevices::hcl.colors(length(config_levels), palette = "Dark 3")
}
config_colours <- setNames(palette[seq_along(config_levels)], config_levels)
configuration_key$colour <- unname(config_colours[config_levels])
write.csv(configuration_key, file.path(output_dir, "configuration_colours.csv"), row.names = FALSE)

set_configuration_order <- function(data) {
  if (any(!data$config_id %in% config_levels)) stop("CSV contains an unknown configuration ID.")
  data$config_id <- factor(data$config_id, levels = config_levels)
  data
}
run_summary <- set_configuration_order(run_summary)
configuration_summary <- set_configuration_order(configuration_summary)
paired <- set_configuration_order(paired)
allocation_results <- set_configuration_order(allocation_results)

colour_scales <- function() list(
  scale_colour_manual(values = config_colours, drop = FALSE, name = "Configuration"),
  scale_fill_manual(values = config_colours, drop = FALSE, name = "Configuration")
)
plot_theme <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "#EDF2F7", colour = "#CFD8E3"),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"), legend.position = "none")

# Hide the redundant screening strip when only one phase exists.
# When multiple phases exist, retain separate panels to avoid pooling results.
phase_facets <- function(data) {
  if (length(unique(data$phase)) > 1) facet_wrap(~phase) else NULL
}
start_facets <- function(data) {
  if (length(unique(data$phase)) > 1) facet_wrap(~phase + starting_suggestion)
  else facet_wrap(~starting_suggestion)
}

# ------------------------------------------------------------
# Variance and paired differences
# ------------------------------------------------------------
variance_plot <- ggplot(filter(run_summary, feasible),
  aes(config_id, variance, fill = config_id, colour = config_id)) +
  geom_boxplot(alpha = 0.35, outlier.shape = NA) +
  geom_point(alpha = 0.65, size = 1.5,
             position = position_jitter(width = 0.12, height = 0, seed = 42)) +
  colour_scales() + phase_facets(run_summary) + coord_flip() + plot_theme +
  labs(title = "Variance by configuration", x = "Configuration",
       y = "Variance (feasible runs only)")
variance_plot

paired_variance_plot <- ggplot(paired,
  aes(config_id, variance_change_pct, fill = config_id, colour = config_id)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey45") +
  geom_boxplot(alpha = 0.35, outlier.shape = NA, na.rm = TRUE) +
  geom_point(alpha = 0.65, size = 1.5, na.rm = TRUE,
             position = position_jitter(width = 0.12, height = 0, seed = 42)) +
  colour_scales() + start_facets(paired) + coord_flip() + plot_theme +
  labs(title = "Paired comparison with baseline", x = "Configuration",
       y = "Variance change from baseline (%)", subtitle = "Negative values favour the configuration; feasible pairs only")
paired_variance_plot

# ------------------------------------------------------------
# Feasibility and fallback rates
# ------------------------------------------------------------
# Dots make zero rates visible. NA labels distinguish unrecorded rates from 0%.
rate_plot <- function(column, title, axis_label) {
  data <- configuration_summary
  data$rate <- data[[column]]
  data$rate_label <- ifelse(is.na(data$rate), "NA", sprintf("%.0f%%", 100 * data$rate))
  data$label_position <- ifelse(is.na(data$rate), 0, data$rate)
  ggplot(data, aes(config_id, rate, fill = config_id, colour = config_id)) +
    geom_col(alpha = 0.8, width = 0.7, na.rm = TRUE) +
    geom_point(size = 2.2, na.rm = TRUE) +
    geom_text(aes(y = label_position, label = rate_label), hjust = -0.25,
              colour = "#303030", size = 3) +
    scale_y_continuous(limits = c(0, 1.12), breaks = seq(0, 1, 0.25),
                        labels = function(x) paste0(100 * x, "%")) +
    colour_scales() + phase_facets(data) + coord_flip() + plot_theme +
    labs(title = title, x = "Configuration", y = axis_label,
         caption = "NA = rate not recorded; it does not mean zero.")
}
feasibility_plot <- rate_plot("feasibility_rate", "Feasibility by configuration", "Feasible runs")
fallback_plot <- rate_plot("fallback_rate", "Fallback to the starting allocation", "Fallback among completed runs")
feasibility_plot
fallback_plot

# ------------------------------------------------------------
# Runtime and generations
# ------------------------------------------------------------
# Use a complete legend instead of overlapping or suppressed point labels.
runtime_plot <- ggplot(configuration_summary,
  aes(median_seconds / 3600, median_variance, colour = config_id)) +
  geom_point(size = 3.5, alpha = 0.85, na.rm = TRUE) +
  scale_colour_manual(values = config_colours, drop = FALSE, name = "Configuration") +
  phase_facets(configuration_summary) + plot_theme +
  theme(legend.position = "right") +
  labs(title = "Runtime and solution quality", x = "Median runtime (hours)",
       y = "Median feasible variance")
runtime_plot

generation_plot <- ggplot(run_summary,
  aes(config_id, generations, fill = config_id, colour = config_id)) +
  geom_boxplot(alpha = 0.35, outlier.shape = NA, na.rm = TRUE) +
  geom_point(alpha = 0.65, size = 1.5, na.rm = TRUE,
             position = position_jitter(width = 0.12, height = 0, seed = 42)) +
  geom_point(data = configuration_summary, aes(y = maxiter),
             shape = 4, size = 3, colour = "#333333") +
  colour_scales() + phase_facets(run_summary) + coord_flip() + plot_theme +
  labs(title = "Generations used", x = "Configuration", y = "Actual generations",
       subtitle = "Dark crosses mark the maximum allowed generations")
generation_plot

# ------------------------------------------------------------
# Allocation stability: separate panels for each stratum and start
# ------------------------------------------------------------
allocation_plot_data <- allocation_results |>
  mutate(across(c(c_h, k_h, n_h), ~ifelse(feasible, .x, NA_real_)))
allocation_plot <- function(column, title, axis_label) {
  facets <- if (length(unique(allocation_plot_data$phase)) > 1) {
    facet_wrap(~phase + starting_suggestion + stratum, scales = "free_y", ncol = 3)
  } else {
    facet_wrap(~starting_suggestion + stratum, scales = "free_y", ncol = 3)
  }
  ggplot(allocation_plot_data,
    aes(config_id, .data[[column]], fill = config_id, colour = config_id)) +
    geom_boxplot(alpha = 0.35, outlier.shape = NA, na.rm = TRUE) +
    geom_point(alpha = 0.6, size = 1, na.rm = TRUE,
               position = position_jitter(width = 0.13, height = 0, seed = 42)) +
    colour_scales() + facets + plot_theme +
    theme(axis.text.x = element_text(angle = 65, hjust = 1, size = 8)) +
    labs(title = title, x = "Configuration", y = axis_label,
         subtitle = "Feasible runs only; each panel has its own vertical scale")
}
ch_stability_plot <- allocation_plot("c_h", "Cluster allocation stability", "Sampled clusters (c_h)")
kh_stability_plot <- allocation_plot("k_h", "Within-cluster allocation stability", "Units per sampled cluster (k_h)")
nh_stability_plot <- allocation_plot("n_h", "Stratum sample-size stability", "Sampled units (n_h)")
ch_stability_plot
kh_stability_plot
nh_stability_plot

# ------------------------------------------------------------
# Save all nine figures as PNG and vector PDF
# ------------------------------------------------------------
plots <- list(variance_boxplots = variance_plot, feasibility_rates = feasibility_plot,
              runtime_vs_variance = runtime_plot, paired_variance_change = paired_variance_plot,
              fallback_rates = fallback_plot, generation_boxplots = generation_plot,
              ch_stability = ch_stability_plot, kh_stability = kh_stability_plot,
              nh_stability = nh_stability_plot)
panel_count <- nrow(distinct(allocation_plot_data, phase, starting_suggestion, stratum))
for (name in names(plots)) {
  is_allocation <- name %in% c("ch_stability", "kh_stability", "nh_stability")
  plot_width <- if (is_allocation) 18 else 13
  plot_height <- if (is_allocation) max(11, 3.8 * ceiling(panel_count / 3)) else 8
  for (extension in c("png", "pdf")) {
    ggsave(file.path(output_dir, paste0(name, ".", extension)), plot = plots[[name]],
           width = plot_width, height = plot_height, dpi = plot_dpi,
           bg = "white", limitsize = FALSE)
  }
}
message("Colour plots saved to: ", normalizePath(output_dir))

