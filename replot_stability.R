# Recreate the existing joint-GA diagnostic plots without rerunning the GA.
library(dplyr)
library(tidyr)
library(ggplot2)

input_dir <- 'joint_ga_stability_plots'
output_dir <- 'joint_ga_stability_plots/reformatted'
plot_dpi <- 200

read_results <- function(name) {
  path <- file.path(input_dir, paste0(name, '.csv'))
  if (!file.exists(path)) stop('Missing results file: ', path)
  read.csv(path, stringsAsFactors = FALSE)
}
run_summary <- read_results('run_summary')
suggestion_metric_comparison <- read_results('suggestion_metric_comparison')
ga_vs_basic_variance_comparison <- read_results('ga_vs_basic_variance_comparison')
basic_benchmark_variances <- read_results('basic_benchmark_variances')
mean_variance <- mean(run_summary$variance, na.rm = TRUE)
best_variance <- min(run_summary$variance, na.rm = TRUE)

allocation_file <- file.path(input_dir, 'allocations.csv')
if (file.exists(allocation_file)) allocations <- read.csv(allocation_file)
have_allocations <- exists('allocations') && is.data.frame(allocations) &&
  all(c('starting_suggestion', 'seed', 'stratum', 'c_h', 'k_h', 'n_h') %in% names(allocations))
if (have_allocations) {
  missing_runs <- anti_join(distinct(run_summary, starting_suggestion, seed),
                           distinct(allocations, starting_suggestion, seed),
                           by = c('starting_suggestion', 'seed'))
  extra_runs <- anti_join(distinct(allocations, starting_suggestion, seed),
                         distinct(run_summary, starting_suggestion, seed),
                         by = c('starting_suggestion', 'seed'))
  if (nrow(missing_runs) || nrow(extra_runs)) stop('Allocation run IDs do not match run_summary.csv')
}

# Starting-strategy colours match the original nh stability figure.
# Keep the theme consistent with replot_hyperparameter_sensitivity.R.
suggestion_names <- sort(unique(run_summary$starting_suggestion))
palette <- c('#4477AA', '#EE6677', '#228833', '#AA3377')
suggestion_colours <- c(
  '1. Proportional to stratum size' = '#F8766D',
  '2. Proportional to between-cluster SD' = '#7CAE00',
  '3. Neyman proportional' = '#00BFC4',
  '4. Equal c_h allocation' = '#C77CFF'
)
if (!all(suggestion_names %in% names(suggestion_colours))) {
  stop('Add named colours for unrecognised starting strategies')
}
plot_theme <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = 'bottom',
        strip.background = element_rect(fill = '#EDF2F7'),
        strip.text = element_text(face = 'bold'))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(data.frame(starting_suggestion = suggestion_names, colour = unname(suggestion_colours[suggestion_names])),
          file.path(output_dir, 'starting_suggestion_colours.csv'), row.names = FALSE)

# Original plot types and statistical summaries are retained below.

 variance_seed_plot <- ggplot(run_summary, aes(x = seed, y = variance, colour = starting_suggestion, shape = feasible)) + 
    geom_line(aes(group = starting_suggestion), 
    linewidth = 0.5, alpha = 0.6) + geom_point(size = 2.8) + scale_shape_manual(values = c(`TRUE` = 16, 
    `FALSE` = 4), labels = c(`TRUE` = "Feasible", `FALSE` = "Infeasible")) + labs(title = "Joint GA stability across random seeds.", x = "Random seed", 
    y = "Estimated variance", colour = "Starting suggestion", shape = "Solution") + theme_minimal(base_size = 12) + 
    theme(legend.position = "bottom", panel.grid.minor = element_blank()) + guides(shape = "none")

variance_distribution_plot <- ggplot(run_summary, aes(x = variance, fill = starting_suggestion)) + geom_histogram(bins = 12, 
    boundary = best_variance, alpha = 0.65, position = "identity") + geom_vline(xintercept = mean(run_summary$variance, 
    na.rm = TRUE), linetype = "dashed", linewidth = 0.8) + geom_vline(xintercept = best_variance, linetype = "dotted", 
    linewidth = 0.8) + labs(title = "Distribution of final GA variances", subtitle = paste0("Mean = ", 
    round(mean(run_summary$variance, na.rm = TRUE), 6), "; SD = ", round(sd(run_summary$variance, na.rm = TRUE), 
        6)), x = "Final variance", y = "Number of runs", fill = "Starting suggestion") + theme_minimal(base_size = 12) + 
    theme(panel.grid.minor = element_blank()) 

variance_boxplot <- ggplot(run_summary, aes(x = reorder(starting_suggestion, variance, FUN = median), 
    y = variance, fill = starting_suggestion)) + geom_boxplot(width = 0.3, outlier.shape = NA) + geom_jitter(aes(shape = feasible), 
    width = 0.08, height = 0, size = 2.4, alpha = 0.75) + scale_shape_manual(values = c(`TRUE` = 16, 
    `FALSE` = 4)) + labs(title = "Variation in the final objective value.", 
    x = NULL, y = "Estimated variance", fill = "Starting suggestion", shape = "Feasible") + coord_flip() + 
    theme_minimal(base_size = 12) + theme(legend.position = "bottom", panel.grid.minor = element_blank()) + guides(shape = "none")

suggestion_metric_comparison_plot <- ggplot(suggestion_metric_comparison, aes(x = reorder(starting_suggestion, 
    value, FUN = median), y = value, fill = starting_suggestion)) + geom_col(show.legend = FALSE) + facet_wrap(~metric, 
    scales = "free_x", ncol = 2) + coord_flip() + labs(title = "Joint GA diagnostic comparison by starting suggestion.", 
    subtitle = "Each metric summarizes repeated runs from the same starting suggestion.", x = NULL, y = "Metric value") + 
    theme_minimal(base_size = 12) + theme(strip.text = element_text(face = "bold"), panel.grid.minor = element_blank()) 

ga_vs_basic_variance_plot <- ggplot(ga_vs_basic_variance_comparison, aes(x = reorder(method, mean_variance), 
    y = mean_variance, fill = comparison_type)) + geom_col() + geom_errorbar(aes(ymin = minimum_variance, 
    ymax = maximum_variance), width = 0.2) + coord_flip() + labs(title = "GA results compared with basic allocation benchmarks.", x = NULL, y = "Estimated variance", 
    fill = "Allocation type") + theme_minimal(base_size = 12) + theme(panel.grid.minor = element_blank()) 

 basic_benchmark_variance_histogram <- ggplot(basic_benchmark_variances, aes(x = reorder(benchmark, variance), 
    y = variance, fill = benchmark)) + geom_col(width = 0.7, show.legend = FALSE) + coord_flip() + labs(title = "Basic benchmark variances", 
    subtitle = "Closed-form non-GA allocations used as comparison benchmarks", x = NULL, y = "Estimated variance") + 
    theme_minimal(base_size = 12) + theme(panel.grid.minor = element_blank()) 



 if (have_allocations) {
allocation_long <- dplyr::mutate(tidyr::pivot_longer(dplyr::select(allocations, starting_suggestion, 
    seed, stratum, c_h, k_h, n_h), cols = c(c_h, k_h, n_h), names_to = "allocation_measure", values_to = "value"), 
    allocation_measure = factor(allocation_measure, levels = c("c_h", "k_h", "n_h"), labels = c("Sampled clusters (c_h)", 
        "Units per cluster (k_h)", "Stratum sample size (n_h)")))
} 

 if (have_allocations) {
combined_allocation_plot <- ggplot(allocation_long, aes(x = stratum, y = value, fill = starting_suggestion)) + 
    geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.15, alpha = 0.35, size = 1) + facet_grid(starting_suggestion ~ 
    allocation_measure, scales = "free_y") + labs(title = "Joint allocation stability across random seeds", 
    subtitle = "Narrow boxes indicate similar allocations across runs", x = "Stratum", y = "Allocated value", 
    fill = "Starting suggestion") + coord_flip() + theme_minimal(base_size = 12) + theme(strip.text = element_text(face = "bold"), 
    panel.grid.minor = element_blank())
} 

repair_diagnostics <- dplyr::mutate(tidyr::pivot_longer(dplyr::select(run_summary, starting_suggestion, 
    seed, sample_size_gap, repair_iterations), cols = c(sample_size_gap, repair_iterations), names_to = "diagnostic", 
    values_to = "value"), diagnostic = factor(diagnostic, levels = c("sample_size_gap", "repair_iterations"), 
    labels = c("Final sample-size gap", "Repair iterations"))) 

repair_diagnostic_plot <- ggplot(repair_diagnostics, aes(x = seed, y = value, colour = starting_suggestion)) + 
    geom_line(aes(group = starting_suggestion), alpha = 0.6) + geom_point(size = 2) + facet_wrap(~diagnostic, 
    scales = "free_y", ncol = 1) + labs(title = "Repair behaviour across random seeds.", x = "Random seed", 
    y = NULL, colour = "Starting suggestion") + theme_minimal(base_size = 12) + theme(strip.text = element_text(face = "bold"), 
    panel.grid.minor = element_blank()) 

# ------------------------------------------------------------
# Apply shared styling and export PNG and PDF copies
# ------------------------------------------------------------
plots <- list('01_variance_by_seed' = variance_seed_plot,
              '02_variance_distribution' = variance_distribution_plot,
              '03_variance_boxplot' = variance_boxplot,
              '04_starting_suggestion_metric_comparison' = suggestion_metric_comparison_plot,
              '04b_ga_vs_basic_variance_comparison' = ga_vs_basic_variance_plot,
              '04c_basic_benchmark_variance_histogram' = basic_benchmark_variance_histogram,
              '10_repair_diagnostics' = repair_diagnostic_plot)
if (have_allocations) {
  plots[['05_ch_stability_by_stratum']] <- c_h_stability_plot
  plots[['06_kh_stability_by_stratum']] <- k_h_stability_plot
  plots[['07_nh_stability_by_stratum']] <- n_h_stability_plot
  plots[['08_combined_allocation_stability']] <- combined_allocation_plot
}
# Preserve the meaning of benchmark colours: these encode comparison types
# or benchmark names, rather than starting strategies.
benchmark_colours <- setNames(palette[seq_along(unique(basic_benchmark_variances$benchmark))],
                              sort(unique(basic_benchmark_variances$benchmark)))
type_colours <- setNames(palette[seq_along(unique(ga_vs_basic_variance_comparison$comparison_type))],
                         sort(unique(ga_vs_basic_variance_comparison$comparison_type)))
for (name in names(plots)) {
  p <- plots[[name]] + plot_theme
  # Translucent boxes keep median lines visible, matching sensitivity plots.
  for (j in seq_along(p$layers)) {
    if (inherits(p$layers[[j]]$geom, 'GeomBoxplot')) {
      p$layers[[j]]$aes_params$alpha <- 0.4
    }
  }
  if (name == '04b_ga_vs_basic_variance_comparison') {
    p <- p + scale_fill_manual(values = type_colours)
  } else if (name == '04c_basic_benchmark_variance_histogram') {
    p <- p + scale_fill_manual(values = benchmark_colours)
  } else {
    if (!is.null(p$mapping$fill)) p <- p + scale_fill_manual(values = suggestion_colours)
    if (!is.null(p$mapping$colour)) p <- p + scale_colour_manual(values = suggestion_colours)
    # Colour point/box outlines by the same strategy as the fill.
    if (!is.null(p$mapping$fill) && name != '04_starting_suggestion_metric_comparison') {
      p <- p + aes(colour = starting_suggestion) +
        scale_colour_manual(values = suggestion_colours, name = 'Starting suggestion')
    }
  }
  # The original combined figure has horizontal numeric axes after coord_flip.
  # Free the numeric scale so c_h, k_h and n_h remain legible separately.
  if (name == '08_combined_allocation_stability') {
    p <- p + facet_grid(starting_suggestion ~ allocation_measure, scales = 'free')
  }
  p <- p + guides(fill = guide_legend(nrow = 2, byrow = TRUE),
                   colour = guide_legend(nrow = 2, byrow = TRUE))
  plots[[name]] <- p
  height <- if (name == '08_combined_allocation_stability') 12 else 8
  for (extension in c('png', 'pdf')) {
    ggsave(file.path(output_dir, paste0(name, '.', extension)), p,
           width = 13, height = height, dpi = plot_dpi, bg = 'white')
  }
}
if (!have_allocations) {
  message('Allocation plots 05-08 were not recreated: per-seed allocations are missing. ',
          'Load the original allocations object or save it to ', allocation_file,
          ', then rerun this script. Summary statistics cannot reconstruct boxplots.')
}
message('Saved ', length(plots), ' reformatted figures as PNG and PDF in ', normalizePath(output_dir))
# For interactive inspection, e.g.: plots[['01_variance_by_seed']]

