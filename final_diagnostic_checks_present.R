rm(list = ls())

library(dplyr)
library(readr)
library(GA)
library(ggplot2)
library(tidyr)

joint_ga_pop_size <- 200
joint_ga_maxiter <- 500
joint_ga_run <- 50
joint_ga_seed_count <- 100 # The number of multiple runs - increase to at least 25 for final check for presentation
# Because there are 4 c_h starting solutions, there will be 4 * joint_ga_seed_count runs of the GA
joint_ga_max_candidate_c <- 1000

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

# Read toy data
toy_data <- read.csv("ecdc-2021-public-v2.csv")

# Convert outcome variable to numeric
toy_data$count_children_present <- readr::parse_number(
  toy_data$count_children_present
)

summarise_joint_ga_result <- function(result, starting_suggestion, seed) {
  data.frame(
    starting_suggestion = starting_suggestion,
    seed = seed,
    variance = result$variance,
    total_n = result$total_n,
    total_c = sum(result$allocation$c_h),
    total_cost = result$total_cost,
    feasible = result$feasible,
    sample_size_gap = result$repair$sample_size_gap,
    repair_iterations = result$repair$iterations
  )
}

add_joint_allocation_context <- function(allocation, starting_suggestion, seed, strata) {
  allocation |>
    dplyr::mutate(
      starting_suggestion = starting_suggestion,
      seed = seed,
      W_h = strata$W_h,
      C_h = strata$C_h,
      N_h = strata$N_h,
      n_share = n_h / sum(n_h),
      population_share = W_h,
      cluster_population_share = C_h / sum(C_h),
      sampled_cluster_share = c_h / sum(c_h),
      cluster_sampling_fraction = c_h / C_h
    )
}

# SUMMARY HELPER SO THAT WE CAN COMPARE THE METHODS OF ALLOCATING C_H
summarise_basic_allocation <- function(benchmark_name,
                                       c_h,
                                       k_h,
                                       strata,
                                       total_n,
                                       max_cost,
                                       lower_c,
                                       upper_c,
                                       lower_k,
                                       upper_k) {
  repaired <- repair_ckh(
    c_h = c_h,
    k_h = k_h,
    total_n = total_n,
    lower_c = lower_c,
    upper_c = upper_c,
    lower_k = lower_k,
    upper_k = upper_k,
    cost_cluster = strata$cost_cluster,
    cost_unit = strata$cost_unit
  )
  
  allocation <- data.frame(
    benchmark = benchmark_name,
    stratum = strata$stratum,
    c_h = repaired$c_h,
    k_h = repaired$k_h,
    n_h = repaired$c_h * repaired$k_h
  )
  
  total_cost <- cost_twostage(
    c_h = allocation$c_h,
    k_h = allocation$k_h,
    cost_cluster = strata$cost_cluster,
    cost_unit = strata$cost_unit
  )
  
  constraint_check <- check_constraints(
    c_h = allocation$c_h,
    k_h = allocation$k_h,
    total_n = total_n,
    max_cost = max_cost,
    cost_cluster = strata$cost_cluster,
    cost_unit = strata$cost_unit,
    min_clusters = 1L,
    min_units = 1L,
    max_clusters = upper_c,
    max_units = upper_k
  )
  
  variance <- variance_twostage(
    c_h = allocation$c_h,
    k_h = allocation$k_h,
    W_h = strata$W_h,
    sigma2_bh = strata$sigma2_bh,
    sigma2_wh = strata$sigma2_wh
  )
  
  summary <- data.frame(
    comparison_type = "Basic allocation",
    method = benchmark_name,
    number_of_runs = 1L,
    feasible_runs = as.integer(constraint_check$feasible),
    feasibility_rate = as.numeric(constraint_check$feasible),
    minimum_variance = variance,
    median_variance = variance,
    mean_variance = variance,
    maximum_variance = variance,
    sd_variance = NA_real_,
    cv_variance_pct = NA_real_,
    mean_total_cost = total_cost,
    mean_sample_size_gap = sum(allocation$n_h) - total_n,
    mean_repair_iterations = repaired$iterations
  )
  
  list(
    summary = summary,
    allocation = allocation,
    repair = repaired,
    constraint_check = constraint_check
  )
}

# THE HARD CODED EQUAL, PROPORTIONAL AND NEYMAN ALLOCATIONS
build_basic_allocation_benchmarks <- function(strata,
                                              total_c,
                                              total_n,
                                              max_cost,
                                              lower_c,
                                              upper_c,
                                              lower_k,
                                              upper_k) {
  equal_c_h <- allocate_ch_by_weight(
    weights = rep(1, nrow(strata)), # equal weights
    total_c = total_c,
    upper = upper_c,
    lower = lower_c
  )
  
  equal_k_h <- rep(
    max(1L, floor(total_n / sum(equal_c_h))),
    nrow(strata)
  )
  
  proportional_c_h <- allocate_ch_by_weight(
    weights = strata$N_h,
    total_c = total_c,
    upper = upper_c,
    lower = lower_c
  )
  
  proportional_k_h <- pmax(
    lower_k,
    round((total_n * strata$W_h) / proportional_c_h)
  )
  
  neyman_sd <- sqrt(strata$sigma2_bh + strata$sigma2_wh)
  neyman_weights <- strata$N_h * neyman_sd
  
  neyman_c_h <- allocate_ch_by_weight(
    weights = neyman_weights,
    total_c = total_c,
    upper = upper_c,
    lower = lower_c
  )
  
  neyman_n_h <- pmax(
    lower_c * lower_k,
    round(total_n * neyman_weights / sum(neyman_weights))
  )
  
  neyman_k_h <- pmax(
    lower_k,
    round(neyman_n_h / neyman_c_h)
  )
  
  equal_benchmark <- summarise_basic_allocation(
    benchmark_name = "Basic equal allocation",
    c_h = equal_c_h,
    k_h = equal_k_h,
    strata = strata,
    total_n = total_n,
    max_cost = max_cost,
    lower_c = lower_c,
    upper_c = upper_c,
    lower_k = lower_k,
    upper_k = upper_k
  )
  
  proportional_benchmark <- summarise_basic_allocation(
    benchmark_name = "Basic stratum proportional allocation",
    c_h = proportional_c_h,
    k_h = proportional_k_h,
    strata = strata,
    total_n = total_n,
    max_cost = max_cost,
    lower_c = lower_c,
    upper_c = upper_c,
    lower_k = lower_k,
    upper_k = upper_k
  )
  
  neyman_benchmark <- summarise_basic_allocation(
    benchmark_name = "Basic Neyman allocation",
    c_h = neyman_c_h,
    k_h = neyman_k_h,
    strata = strata,
    total_n = total_n,
    max_cost = max_cost,
    lower_c = lower_c,
    upper_c = upper_c,
    lower_k = lower_k,
    upper_k = upper_k
  )
  
  list(
    summary = dplyr::bind_rows(
      equal_benchmark$summary,
      proportional_benchmark$summary,
      neyman_benchmark$summary
    ),
    allocations = dplyr::bind_rows(
      equal_benchmark$allocation,
      proportional_benchmark$allocation,
      neyman_benchmark$allocation
    )
  )
}

# ------------------------------------------------------------
# Full ECD data joint-allocation test across all provinces
# ------------------------------------------------------------

# Prepare full strata-level data using the sampling frame workflow
survey_data <- read_survey_data(
  data = toy_data,
  stratum_id = "province",
  cluster_id = "id_ward",
  y_var = "count_children_present",
  na_action = "warn"
)

total_n <- 5000
max_cost <- 1400000
unit_cost <- 100
cluster_cost <- 2500

frame_c_h <- frame_c(survey_data, max_cost, unit_cost, cluster_cost, total_n)

# Produces a proportional allocation of the available clusters, bounded by cost
full_frame <- prepare_sampling_frame(
  data = survey_data,
  y_var = "count_children_present",
  total_c = frame_c_h,
  cost_cluster = cluster_cost,
  cost_unit = unit_cost
)

strata <- full_frame$strata

strata # check

L <- nrow(strata)

lower_c <- rep(1L, L)
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

total_c <- joint_start$total_c
total_c 

# Build the four starting suggestions used in the diagnostic comparisons:
# three proportional c_h allocations and one equal c_h allocation. These are
# starting suggestions only; the joint GA can still move both c_h and k_h.
joint_starting_suggestions <- joint_start$all_starts |>
  dplyr::rename(
    repaired_sample_size_gap = sample_size_gap,
    repaired_feasible = feasible
  )

starting_suggestion_diagnostics <- joint_starting_suggestions |>
  dplyr::group_by(starting_suggestion) |>
  dplyr::summarise(
    total_initial_c = sum(initial_c_h),
    total_repaired_c = sum(repaired_c_h),
    total_repaired_n = dplyr::first(repaired_total_n),
    total_repaired_cost = dplyr::first(repaired_total_cost),
    sample_size_gap = dplyr::first(repaired_sample_size_gap),
    repaired_exactly = dplyr::first(repaired_exactly),
    repaired_feasible = dplyr::first(repaired_feasible),
    repair_iterations = dplyr::first(repair_iterations),
    .groups = "drop"
  )

#joint_starting_suggestions
#starting_suggestion_diagnostics

# default_start <- joint_starting_suggestions |>
#   dplyr::filter(starting_suggestion == "3. Neyman proportional")

# One run with Neyman proportional c_h allocation 
#ga_result <- allocate_twostage_ga_joint(
#  strata = strata,
#  total_n = total_n,
#  max_cost = max_cost,
#  popSize = joint_ga_pop_size,
#  maxiter = joint_ga_maxiter,
#  run = joint_ga_run,
#  seed = 123,
#  suggestion = c(default_start$repaired_c_h, default_start$repaired_k_h)
#)

# ga_result$allocation
# ga_result$variance
# ga_result$total_n
# ga_result$total_cost
# ga_result$feasible
# ga_result$constraint_check
# ga_result$repair

################################################################################
#                    MULTIPLE RUN TESTS
################################################################################

set.seed(123)

strata <- strata
total_n <- total_n
max_cost <- max_cost

seeds <- seq_len(joint_ga_seed_count) # Increase/decrease as needed; many seeds can take long

starting_suggestion_names <- unique(joint_starting_suggestions$starting_suggestion)

ga_run_grid <- expand.grid(
  starting_suggestion = starting_suggestion_names,
  seed = seeds,
  stringsAsFactors = FALSE
)

ga_runs <- lapply(seq_len(nrow(ga_run_grid)), function(i) {
  start_i <- joint_starting_suggestions |>
    dplyr::filter(starting_suggestion == ga_run_grid$starting_suggestion[i])
  
  allocate_twostage_ga_joint(
    strata = strata,
    total_n = total_n,
    max_cost = max_cost,
    popSize = joint_ga_pop_size,
    maxiter = joint_ga_maxiter,
    run = joint_ga_run,
    seed = ga_run_grid$seed[i],
    suggestion = c(start_i$repaired_c_h, start_i$repaired_k_h)
  )
})

run_summary <- dplyr::bind_rows(
  lapply(seq_along(ga_runs), function(i) {
    summarise_joint_ga_result(
      result = ga_runs[[i]],
      starting_suggestion = ga_run_grid$starting_suggestion[i],
      seed = ga_run_grid$seed[i]
    )
  })
)

run_summary

#summary(run_summary$variance)

#min(run_summary$variance)
#mean(run_summary$variance)
#sd(run_summary$variance) # should be small if the solution is stable
#sum(run_summary$feasible)

suggestion_run_summary <- run_summary |>
  dplyr::group_by(starting_suggestion) |>
  dplyr::summarise(
    number_of_runs = dplyr::n(),
    feasible_runs = sum(feasible, na.rm = TRUE),
    feasibility_rate = mean(feasible, na.rm = TRUE),
    minimum_variance = min(variance, na.rm = TRUE),
    median_variance = median(variance, na.rm = TRUE),
    mean_variance = mean(variance, na.rm = TRUE),
    maximum_variance = max(variance, na.rm = TRUE),
    sd_variance = sd(variance, na.rm = TRUE),
    cv_variance_pct =
      100 * sd(variance, na.rm = TRUE) /
      mean(variance, na.rm = TRUE),
    mean_total_cost = mean(total_cost, na.rm = TRUE),
    mean_sample_size_gap = mean(sample_size_gap, na.rm = TRUE),
    mean_repair_iterations = mean(repair_iterations, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::arrange(minimum_variance)

suggestion_run_summary

basic_allocation_benchmarks <- build_basic_allocation_benchmarks(
  strata = strata,
  total_c = total_c,
  total_n = total_n,
  max_cost = max_cost,
  lower_c = lower_c,
  upper_c = upper_c,
  lower_k = lower_k,
  upper_k = upper_k
)

basic_allocation_summary <- basic_allocation_benchmarks$summary
basic_allocation_allocations <- basic_allocation_benchmarks$allocations

basic_benchmark_variances <- basic_allocation_summary |>
  dplyr::transmute(
    benchmark = method,
    variance = mean_variance,
    total_cost = mean_total_cost,
    sample_size_gap = mean_sample_size_gap,
    feasible = feasible_runs == 1L
  ) |>
  dplyr::arrange(variance)

ga_vs_basic_variance_comparison <- dplyr::bind_rows(
  suggestion_run_summary |>
    dplyr::transmute(
      comparison_type = "GA starting suggestion",
      method = starting_suggestion,
      number_of_runs = number_of_runs,
      feasible_runs = feasible_runs,
      feasibility_rate = feasibility_rate,
      minimum_variance = minimum_variance,
      median_variance = median_variance,
      mean_variance = mean_variance,
      maximum_variance = maximum_variance,
      sd_variance = sd_variance,
      cv_variance_pct = cv_variance_pct,
      mean_total_cost = mean_total_cost,
      mean_sample_size_gap = mean_sample_size_gap,
      mean_repair_iterations = mean_repair_iterations
    ),
  basic_allocation_summary
) |>
  dplyr::mutate(
    improvement_vs_basic_equal_pct =
      100 *
      (
        mean_variance[method == "Basic equal allocation"][1] -
          mean_variance
      ) /
      mean_variance[method == "Basic equal allocation"][1],
    improvement_vs_basic_proportional_pct =
      100 *
      (
        mean_variance[method == "Basic stratum proportional allocation"][1] -
          mean_variance
      ) /
      mean_variance[method == "Basic stratum proportional allocation"][1],
    improvement_vs_basic_neyman_pct =
      100 *
      (
        mean_variance[method == "Basic Neyman allocation"][1] -
          mean_variance
      ) /
      mean_variance[method == "Basic Neyman allocation"][1]
  ) |>
  dplyr::arrange(mean_variance)

#basic_allocation_summary
#basic_allocation_allocations
#basic_benchmark_variances
#ga_vs_basic_variance_comparison

# Best GA run
best_index <- which.min(run_summary$variance)

best_run <- ga_runs[[best_index]]

best_run$allocation
best_run$variance
best_run$total_n
best_run$total_cost
best_run$feasible

# Average allocations across runs
allocations <- dplyr::bind_rows(
  lapply(seq_along(ga_runs), function(i) {
    add_joint_allocation_context(
      allocation = ga_runs[[i]]$allocation,
      starting_suggestion = ga_run_grid$starting_suggestion[i],
      seed = ga_run_grid$seed[i],
      strata = strata
    )
  })
)

allocation_summary <- allocations |>
  dplyr::group_by(starting_suggestion, stratum) |>
  dplyr::summarise(
    mean_c_h = mean(c_h),
    sd_c_h = sd(c_h),
    min_c_h = min(c_h),
    max_c_h = max(c_h),
    mean_k_h = mean(k_h),
    sd_k_h = sd(k_h),
    min_k_h = min(k_h),
    max_k_h = max(k_h),
    mean_n_h = mean(n_h),
    sd_n_h = sd(n_h),
    min_n_h = min(n_h),
    max_n_h = max(n_h),
    .groups = "drop"
  )

#allocation_summary

################################################################################
#                         STABILITY PLOTS
################################################################################

stability_output_dir <- "joint_ga_stability_plots"
dir.create(
  stability_output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

# Add useful stability measures
best_variance <- min(run_summary$variance, na.rm = TRUE)
mean_variance <- mean(run_summary$variance, na.rm = TRUE)

run_summary <- run_summary |>
  dplyr::group_by(starting_suggestion) |>
  dplyr::mutate(
    best_suggestion_variance = min(variance, na.rm = TRUE),
    suggestion_variance_gap = variance - best_suggestion_variance,
    suggestion_relative_gap_pct =
      100 * (variance - best_suggestion_variance) / best_suggestion_variance,
    global_variance_gap = variance - best_variance,
    global_relative_gap_pct = 100 * (variance - best_variance) / best_variance,
    run = dplyr::row_number()
  ) |>
  dplyr::ungroup()

# Numerical stability summary
stability_summary <- run_summary |>
  dplyr::group_by(starting_suggestion) |>
  dplyr::summarise(
    number_of_runs = dplyr::n(),
    feasible_runs = sum(feasible, na.rm = TRUE),
    feasibility_rate = mean(feasible, na.rm = TRUE),
    minimum_variance = min(variance, na.rm = TRUE),
    median_variance = median(variance, na.rm = TRUE),
    mean_variance = mean(variance, na.rm = TRUE),
    maximum_variance = max(variance, na.rm = TRUE),
    sd_variance = sd(variance, na.rm = TRUE),
    cv_variance_pct =
      100 * sd(variance, na.rm = TRUE) /
      mean(variance, na.rm = TRUE),
    maximum_relative_gap_pct =
      max(suggestion_relative_gap_pct, na.rm = TRUE),
    best_global_relative_gap_pct =
      min(global_relative_gap_pct, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::arrange(minimum_variance)

stability_summary

suggestion_metric_comparison <- stability_summary |>
  dplyr::select(
    starting_suggestion,
    feasibility_rate,
    minimum_variance,
    median_variance,
    mean_variance,
    sd_variance,
    cv_variance_pct,
    maximum_relative_gap_pct
  ) |>
  tidyr::pivot_longer(
    cols = -starting_suggestion,
    names_to = "metric",
    values_to = "value"
  ) |>
  dplyr::mutate(
    metric = factor(
      metric,
      levels = c(
        "minimum_variance",
        "median_variance",
        "mean_variance",
        "sd_variance",
        "cv_variance_pct",
        "maximum_relative_gap_pct",
        "feasibility_rate"
      ),
      labels = c(
        "Minimum variance",
        "Median variance",
        "Mean variance",
        "SD of variance",
        "CV of variance (%)",
        "Maximum relative gap (%)",
        "Feasibility rate"
      )
    )
  )

suggestion_metric_comparison


# ============================================================
# 1. Variance obtained for every seed
# ============================================================

variance_seed_plot <- ggplot(
  run_summary,
  aes(
    x = seed,
    y = variance,
    colour = starting_suggestion,
    shape = feasible
  )
) +
  geom_hline(
    yintercept = mean_variance,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  geom_line(
    aes(group = starting_suggestion),
    linewidth = 0.5,
    alpha = 0.6
  ) +
  geom_point(size = 2.8) +
  scale_shape_manual(
    values = c(
      "TRUE" = 16,
      "FALSE" = 4
    ),
    labels = c(
      "TRUE" = "Feasible",
      "FALSE" = "Infeasible"
    )
  ) +
  labs(
    title = "Joint GA stability across random seeds",
    subtitle = paste0(
      "Dashed line shows the mean variance: ",
      round(mean_variance, 6)
    ),
    x = "Random seed",
    y = "Estimated variance",
    colour = "Starting suggestion",
    shape = "Solution"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

variance_seed_plot

ggsave(
  filename = file.path(
    stability_output_dir,
    "01_variance_by_seed.png"
  ),
  plot = variance_seed_plot,
  width = 9,
  height = 5.5,
  dpi = 300
)


# ============================================================
# 2. Distribution of the final variance across runs
# ============================================================

variance_distribution_plot <- ggplot(
  run_summary,
  aes(
    x = variance,
    fill = starting_suggestion
  )
) +
  geom_histogram(
    bins = 12,
    boundary = best_variance,
    alpha = 0.65,
    position = "identity"
  ) +
  geom_vline(
    xintercept = mean(run_summary$variance, na.rm = TRUE),
    linetype = "dashed",
    linewidth = 0.8
  ) +
  geom_vline(
    xintercept = best_variance,
    linetype = "dotted",
    linewidth = 0.8
  ) +
  labs(
    title = "Distribution of final GA variances",
    subtitle = paste0(
      "Mean = ",
      round(mean(run_summary$variance, na.rm = TRUE), 6),
      "; SD = ",
      round(sd(run_summary$variance, na.rm = TRUE), 6)
    ),
    x = "Final variance",
    y = "Number of runs",
    fill = "Starting suggestion"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

variance_distribution_plot

ggsave(
  filename = file.path(
    stability_output_dir,
    "02_variance_distribution.png"
  ),
  plot = variance_distribution_plot,
  width = 8,
  height = 5.5,
  dpi = 300
)


# ============================================================
# 3. Boxplot and individual final variance values
# ============================================================

variance_boxplot <- ggplot(
  run_summary,
  aes(
    x = reorder(starting_suggestion, variance, FUN = median),
    y = variance,
    fill = starting_suggestion
  )
) +
  geom_boxplot(
    width = 0.3,
    outlier.shape = NA
  ) +
  geom_jitter(
    aes(shape = feasible),
    width = 0.08,
    height = 0,
    size = 2.4,
    alpha = 0.75
  ) +
  scale_shape_manual(
    values = c(
      "TRUE" = 16,
      "FALSE" = 4
    )
  ) +
  labs(
    title = "Variation in the final objective value",
    subtitle = "Each box summarizes repeated runs from the same starting suggestion",
    x = NULL,
    y = "Estimated variance",
    fill = "Starting suggestion",
    shape = "Feasible"
  ) +
  coord_flip() +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

variance_boxplot

ggsave(
  filename = file.path(
    stability_output_dir,
    "03_variance_boxplot.png"
  ),
  plot = variance_boxplot,
  width = 6.5,
  height = 5.5,
  dpi = 300
)


# ============================================================
# 4. Stability metrics compared across starting suggestions
# ============================================================

suggestion_metric_comparison_plot <- ggplot(
  suggestion_metric_comparison,
  aes(
    x = reorder(starting_suggestion, value, FUN = median),
    y = value,
    fill = starting_suggestion
  )
) +
  geom_col(show.legend = FALSE) +
  facet_wrap(
    ~ metric,
    scales = "free_x",
    ncol = 2
  ) +
  coord_flip() +
  labs(
    title = "Joint GA diagnostic comparison by starting suggestion",
    subtitle = "Each metric summarizes repeated runs from the same starting suggestion",
    x = NULL,
    y = "Metric value"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

suggestion_metric_comparison_plot

ggsave(
  filename = file.path(
    stability_output_dir,
    "04_starting_suggestion_metric_comparison.png"
  ),
  plot = suggestion_metric_comparison_plot,
  width = 11,
  height = 8,
  dpi = 300
)


# ============================================================
# 4b. GA results compared with basic allocation benchmarks
# ============================================================

ga_vs_basic_variance_plot <- ggplot(
  ga_vs_basic_variance_comparison,
  aes(
    x = reorder(method, mean_variance),
    y = mean_variance,
    fill = comparison_type
  )
) +
  geom_col() +
  geom_errorbar(
    aes(
      ymin = minimum_variance,
      ymax = maximum_variance
    ),
    width = 0.2
  ) +
  coord_flip() +
  labs(
    title = "GA results compared with basic allocation benchmarks",
    subtitle = "Bars show mean variance; error bars show min and max across GA runs",
    x = NULL,
    y = "Estimated variance",
    fill = "Allocation type"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ga_vs_basic_variance_plot

ggsave(
  filename = file.path(
    stability_output_dir,
    "04b_ga_vs_basic_variance_comparison.png"
  ),
  plot = ga_vs_basic_variance_plot,
  width = 10,
  height = 6,
  dpi = 300
)


# ============================================================
# 4c. Basic benchmark variance histogram
# ============================================================

basic_benchmark_variance_histogram <- ggplot(
  basic_benchmark_variances,
  aes(
    x = reorder(benchmark, variance),
    y = variance,
    fill = benchmark
  )
) +
  geom_col(width = 0.7, show.legend = FALSE) +
  coord_flip() +
  labs(
    title = "Basic benchmark variances",
    subtitle = "Closed-form non-GA allocations used as comparison benchmarks",
    x = NULL,
    y = "Estimated variance"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

basic_benchmark_variance_histogram

ggsave(
  filename = file.path(
    stability_output_dir,
    "04c_basic_benchmark_variance_histogram.png"
  ),
  plot = basic_benchmark_variance_histogram,
  width = 9,
  height = 5,
  dpi = 300
)


# ============================================================
# 5. Distribution of c_h across runs for each stratum
# ============================================================

c_h_stability_plot <- ggplot(
  allocations,
  aes(
    x = reorder(stratum, c_h, FUN = median),
    y = c_h,
    fill = starting_suggestion
  )
) +
  geom_boxplot(
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.15,
    alpha = 0.45,
    size = 1.5
  ) +
  labs(
    title = "Stability of sampled-cluster allocations",
    subtitle = "Each point is the allocation from one random seed",
    x = "Stratum",
    y = expression(c[h]),
    fill = "Starting suggestion"
  ) +
  facet_wrap(~ starting_suggestion) +
  coord_flip() +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

c_h_stability_plot

ggsave(
  filename = file.path(
    stability_output_dir,
    "05_ch_stability_by_stratum.png"
  ),
  plot = c_h_stability_plot,
  width = 9,
  height = 6,
  dpi = 300
)


# ============================================================
# 6. Distribution of k_h across runs for each stratum
# ============================================================

k_h_stability_plot <- ggplot(
  allocations,
  aes(
    x = reorder(stratum, k_h, FUN = median),
    y = k_h,
    fill = starting_suggestion
  )
) +
  geom_boxplot(
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.15,
    alpha = 0.45,
    size = 1.5
  ) +
  labs(
    title = "Stability of within-cluster sample sizes",
    subtitle = "Each point is the allocation from one random seed",
    x = "Stratum",
    y = expression(k[h]),
    fill = "Starting suggestion"
  ) +
  facet_wrap(~ starting_suggestion) +
  coord_flip() +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

k_h_stability_plot

ggsave(
  filename = file.path(
    stability_output_dir,
    "06_kh_stability_by_stratum.png"
  ),
  plot = k_h_stability_plot,
  width = 9,
  height = 6,
  dpi = 300
)

# ============================================================
# 7. Distribution of n_h = c_h * k_h across runs
# ============================================================

n_h_stability_plot <- ggplot(
  allocations,
  aes(
    x = reorder(stratum, n_h, FUN = median),
    y = n_h,
    fill = starting_suggestion
  )
) +
  geom_boxplot(
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.15,
    alpha = 0.45,
    size = 1.5
  ) +
  labs(
    title = "Stability of total stratum sample sizes",
    subtitle = expression(
      paste(
        "Each point is one GA run;  ",
        n[h] == c[h] * k[h]
      )
    ),
    x = "Stratum",
    y = expression(n[h]),
    fill = "Starting suggestion"
  ) +
  facet_wrap(~ starting_suggestion) +
  coord_flip() +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

n_h_stability_plot

ggsave(
  filename = file.path(
    stability_output_dir,
    "07_nh_stability_by_stratum.png"
  ),
  plot = n_h_stability_plot,
  width = 9,
  height = 6,
  dpi = 300
)


# ============================================================
# 8. Combined faceted allocation stability plot
# ============================================================

allocation_long <- allocations |>
  dplyr::select(
    starting_suggestion,
    seed,
    stratum,
    c_h,
    k_h,
    n_h
  ) |>
  tidyr::pivot_longer(
    cols = c(c_h, k_h, n_h),
    names_to = "allocation_measure",
    values_to = "value"
  ) |>
  dplyr::mutate(
    allocation_measure = factor(
      allocation_measure,
      levels = c("c_h", "k_h", "n_h"),
      labels = c(
        "Sampled clusters (c_h)",
        "Units per cluster (k_h)",
        "Stratum sample size (n_h)"
      )
    )
  )

combined_allocation_plot <- ggplot(
  allocation_long,
  aes(
    x = stratum,
    y = value,
    fill = starting_suggestion
  )
) +
  geom_boxplot(
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.15,
    alpha = 0.35,
    size = 1
  ) +
  facet_grid(
    starting_suggestion ~ allocation_measure,
    scales = "free_y"
  ) +
  labs(
    title = "Joint allocation stability across random seeds",
    subtitle = "Narrow boxes indicate similar allocations across runs",
    x = "Stratum",
    y = "Allocated value",
    fill = "Starting suggestion"
  ) +
  coord_flip() +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

combined_allocation_plot

ggsave(
  filename = file.path(
    stability_output_dir,
    "08_combined_allocation_stability.png"
  ),
  plot = combined_allocation_plot,
  width = 10,
  height = 12,
  dpi = 300
)

# ============================================================
# 9. Feasibility and repair diagnostics across seeds
# ============================================================

repair_diagnostics <- run_summary |>
  dplyr::select(
    starting_suggestion,
    seed,
    sample_size_gap,
    repair_iterations
  ) |>
  tidyr::pivot_longer(
    cols = c(sample_size_gap, repair_iterations),
    names_to = "diagnostic",
    values_to = "value"
  ) |>
  dplyr::mutate(
    diagnostic = factor(
      diagnostic,
      levels = c(
        "sample_size_gap",
        "repair_iterations"
      ),
      labels = c(
        "Final sample-size gap",
        "Repair iterations"
      )
    )
  )

repair_diagnostic_plot <- ggplot(
  repair_diagnostics,
  aes(
    x = seed,
    y = value,
    colour = starting_suggestion
  )
) +
  geom_line(aes(group = starting_suggestion), alpha = 0.6) +
  geom_point(size = 2) +
  facet_wrap(
    ~ diagnostic,
    scales = "free_y",
    ncol = 1
  ) +
  labs(
    title = "Repair behaviour across random seeds",
    x = "Random seed",
    y = "Diagnostic value",
    colour = "Starting suggestion"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

repair_diagnostic_plot

ggsave(
  filename = file.path(
    stability_output_dir,
    "10_repair_diagnostics.png"
  ),
  plot = repair_diagnostic_plot,
  width = 9,
  height = 7,
  dpi = 300
)


# ============================================================
# Save the numerical results
# ============================================================

write.csv(
  starting_suggestion_diagnostics,
  file.path(
    stability_output_dir,
    "starting_suggestion_diagnostics.csv"
  ),
  row.names = FALSE
)

write.csv(
  run_summary,
  file.path(
    stability_output_dir,
    "run_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  allocation_summary,
  file.path(
    stability_output_dir,
    "allocation_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  suggestion_run_summary,
  file.path(
    stability_output_dir,
    "suggestion_run_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  basic_allocation_summary,
  file.path(
    stability_output_dir,
    "basic_allocation_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  basic_allocation_allocations,
  file.path(
    stability_output_dir,
    "basic_allocation_allocations.csv"
  ),
  row.names = FALSE
)

write.csv(
  basic_benchmark_variances,
  file.path(
    stability_output_dir,
    "basic_benchmark_variances.csv"
  ),
  row.names = FALSE
)

write.csv(
  ga_vs_basic_variance_comparison,
  file.path(
    stability_output_dir,
    "ga_vs_basic_variance_comparison.csv"
  ),
  row.names = FALSE
)

write.csv(
  suggestion_metric_comparison,
  file.path(
    stability_output_dir,
    "suggestion_metric_comparison.csv"
  ),
  row.names = FALSE
)

write.csv(
  stability_summary,
  file.path(
    stability_output_dir,
    "stability_summary.csv"
  ),
  row.names = FALSE
)

message(
  "Stability plots and summaries saved to: ",
  stability_output_dir
)

