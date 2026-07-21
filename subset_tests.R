############### TEST WITH TOY DATA ##################################

rm(list = ls())

library(dplyr)
library(readr)
library(GA)

source("variance_twostage.R")
source("cost_twostage.R")
source("check_constraints.R")
source("repair_joint.R")
source("fitness_ch_kh.R")
source("allocate_twostage_joint.R")
source("read_survey_data.R")
source("sampling_frame.R")

# Read toy data
toy_data <- read.csv("ecdc-2021-public-v2.csv")

# Convert outcome variable to numeric
toy_data$count_children_present <- readr::parse_number(
  toy_data$count_children_present
)

allocate_ch_by_weight <- function(weights,
                                  total_c,
                                  upper,
                                  lower = rep(1L, length(weights))) {
  weights <- as.numeric(weights)
  upper <- as.integer(upper)
  lower <- as.integer(lower)
  total_c <- as.integer(total_c)
  
  if (length(lower) != length(weights) || length(upper) != length(weights)) {
    stop("weights, lower, and upper must have the same length.")
  }
  
  if (total_c < sum(lower) || total_c > sum(upper)) {
    stop("total_c is infeasible for these c_h bounds.")
  }
  
  weights[!is.finite(weights) | weights < 0] <- 0
  
  if (sum(weights) == 0) {
    weights <- rep(1, length(weights))
  }
  
  c_h <- lower
  remaining <- total_c - sum(c_h)
  capacity <- upper - c_h
  
  if (remaining > 0L) {
    raw_extra <- remaining * weights / sum(weights)
    extra <- pmin(floor(raw_extra), capacity)
    c_h <- c_h + extra
    
    while (sum(c_h) < total_c) {
      eligible <- which(c_h < upper)
      frac <- raw_extra - floor(raw_extra)
      h <- eligible[order(frac[eligible], weights[eligible], decreasing = TRUE)][1]
      c_h[h] <- c_h[h] + 1L
    }
  }
  
  as.integer(c_h)
}

################################################################################
#                SUBSET TEST
################################################################################

# This reduces the GA decision vector from 2 * 9 provinces to 2 * 3 provinces.
test_provinces <- c("Eastern Cape", "Gauteng", "KwaZulu-Natal")

toy_data_subset <- toy_data |>
  dplyr::filter(province %in% test_provinces)

# Prepare strata-level data using the sampling frame workflow
subset_survey_data <- read_survey_data(
  data = toy_data_subset,
  stratum_id = "province",
  cluster_id = "id_ward",
  y_var = "count_children_present",
  na_action = "warn"
)

# Max affordable number of clusters
subset_total_c <- floor((8000 - 10 * 300) / 100)

subset_frame <- prepare_sampling_frame(
  data = subset_survey_data,
  y_var = "count_children_present",
  total_c = subset_total_c,
  cost_cluster = 100,
  cost_unit = 10
)

strata <- subset_frame$strata

# View prepared strata
strata

# Allocation constraints
total_n <- 300
max_cost <- 8000

# Number of strata
L <- nrow(strata)

# Bounds for c_h and k_h
lower_c <- rep(1L, L)
upper_c <- strata$C_h

lower_k <- rep(1L, L)
upper_k <- rep(total_n, L)

# Start from a variance-proportional c_h suggestion.
initial_c_h <- allocate_ch_by_weight(
  weights = strata$sigma2_bh,
  total_c = subset_total_c,
  upper = upper_c,
  lower = lower_c
)
initial_k_h <- rep(8L, L)

repaired <- repair_ckh(
  c_h = initial_c_h,
  k_h = initial_k_h,
  total_n = total_n,
  lower_c = lower_c,
  upper_c = upper_c,
  lower_k = lower_k,
  upper_k = upper_k
)

repaired
sum(repaired$c_h * repaired$k_h)

# Run the joint GA on the reduced province subset.
ga_result <- allocate_twostage_ga_joint(
  strata = strata,
  total_n = total_n,
  max_cost = max_cost,
  popSize = 100,
  maxiter = 200,
  run = 50,
  seed = 123,
  suggestion = c(repaired$c_h, repaired$k_h)
)

ga_result$allocation
ga_result$variance
ga_result$total_n
ga_result$total_cost
ga_result$feasible
ga_result$constraint_check

# ------------------------------------------------------------
# Full ECD data joint-allocation test across all provinces
# ------------------------------------------------------------

# Prepare full strata-level data using the sampling frame workflow
full_survey_data <- read_survey_data(
  data = toy_data,
  stratum_id = "province",
  cluster_id = "id_ward",
  y_var = "count_children_present",
  na_action = "warn"
)

full_total_n <- 900
full_max_cost <- 18000
full_total_c <- floor((full_max_cost - 10 * full_total_n) / 100)

# Produces a proportional allocation of the 90 clusters
full_frame <- prepare_sampling_frame(
  data = full_survey_data,
  y_var = "count_children_present",
  total_c = full_total_c,
  cost_cluster = 100,
  cost_unit = 10
)

full_strata <- full_frame$strata

full_strata

full_L <- nrow(full_strata)

full_lower_c <- rep(1L, full_L)
full_upper_c <- full_strata$C_h

full_lower_k <- rep(1L, full_L)
full_upper_k <- rep(full_total_n, full_L)

# Start from a variance-proportional c_h allocation and deliberately off-target k_h.
# In the joint GA this is only a starting suggestion; c_h can still move.
full_initial_c_h <- allocate_ch_by_weight(
  weights = full_strata$sigma2_bh,
  total_c = full_total_c,
  upper = full_upper_c,
  lower = full_lower_c
)
full_initial_c_h <- pmin(full_initial_c_h, full_upper_c)
full_initial_k_h <- rep(8L, full_L)

full_repaired <- repair_ckh(
  c_h = full_initial_c_h,
  k_h = full_initial_k_h,
  total_n = full_total_n,
  lower_c = full_lower_c,
  upper_c = full_upper_c,
  lower_k = full_lower_k,
  upper_k = full_upper_k
)

full_repaired
sum(full_repaired$c_h * full_repaired$k_h)

full_ga_result <- allocate_twostage_ga_joint(
  strata = full_strata,
  total_n = full_total_n,
  max_cost = full_max_cost,
  popSize = 100,
  maxiter = 200,
  run = 50,
  seed = 123,
  suggestion = c(full_repaired$c_h, full_repaired$k_h)
)

full_ga_result$allocation
full_ga_result$variance
full_ga_result$total_n
full_ga_result$total_cost
full_ga_result$feasible
full_ga_result$constraint_check
full_ga_result$repair

################################################################################
#                    MULTIPLE RUN TESTS
################################################################################

set.seed(2026)

strata <- full_strata
total_n <- full_total_n
max_cost <- full_max_cost

seeds <- 1:10

starting_c_h <- full_repaired$c_h
starting_k_h <- full_repaired$k_h

starting_c_h
starting_k_h

ga_runs <- lapply(seeds, function(s) {
  allocate_twostage_ga_joint(
    strata = strata,
    total_n = total_n,
    max_cost = max_cost,
    popSize = 100,
    maxiter = 200,
    run = 50,
    seed = s,
    suggestion = c(starting_c_h, starting_k_h)
  )
})

run_summary <- data.frame(
  seed = seeds,
  variance = sapply(ga_runs, function(x) x$variance),
  total_n = sapply(ga_runs, function(x) x$total_n),
  total_cost = sapply(ga_runs, function(x) x$total_cost),
  feasible = sapply(ga_runs, function(x) x$feasible),
  sample_size_gap = sapply(ga_runs, function(x) x$repair$sample_size_gap),
  repair_iterations = sapply(ga_runs, function(x) x$repair$iterations)
)

run_summary

summary(run_summary$variance)

min(run_summary$variance)
mean(run_summary$variance)
sd(run_summary$variance) # should be small if the solution is stable
sum(run_summary$feasible)

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
    ga_runs[[i]]$allocation |>
      dplyr::mutate(seed = seeds[i])
  })
)

allocation_summary <- allocations |>
  dplyr::group_by(stratum) |>
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

allocation_summary

library(ggplot2)
library(tidyr)

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

run_summary <- run_summary |>
  dplyr::mutate(
    variance_gap = variance - best_variance,
    relative_gap_pct = 100 * (variance - best_variance) / best_variance,
    run = dplyr::row_number()
  )

# Numerical stability summary
stability_summary <- run_summary |>
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
      max(relative_gap_pct, na.rm = TRUE)
  )

stability_summary


# ============================================================
# 1. Variance obtained for every seed
# ============================================================

variance_seed_plot <- ggplot(
  run_summary,
  aes(
    x = seed,
    y = variance,
    shape = feasible
  )
) +
  geom_hline(
    yintercept = best_variance,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  geom_line(
    aes(group = 1),
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
      "Dashed line shows the best variance: ",
      round(best_variance, 6)
    ),
    x = "Random seed",
    y = "Estimated variance",
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
  aes(x = variance)
) +
  geom_histogram(
    bins = 12,
    boundary = best_variance
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
    y = "Number of runs"
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
    x = "Joint GA",
    y = variance
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
    subtitle = paste0(
      "Coefficient of variation = ",
      round(stability_summary$cv_variance_pct, 4),
      "%"
    ),
    x = NULL,
    y = "Estimated variance",
    shape = "Feasible"
  ) +
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
# 4. Distribution of c_h across runs for each stratum
# ============================================================

c_h_stability_plot <- ggplot(
  allocations,
  aes(
    x = reorder(stratum, c_h, FUN = median),
    y = c_h
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
    y = expression(c[h])
  ) +
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
# 5. Distribution of k_h across runs for each stratum
# ============================================================

k_h_stability_plot <- ggplot(
  allocations,
  aes(
    x = reorder(stratum, k_h, FUN = median),
    y = k_h
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
    y = expression(k[h])
  ) +
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
# 6. Distribution of n_h = c_h * k_h across runs
# ============================================================

n_h_stability_plot <- ggplot(
  allocations,
  aes(
    x = reorder(stratum, n_h, FUN = median),
    y = n_h
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
    y = expression(n[h])
  ) +
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
# 7. Combined faceted allocation stability plot
# ============================================================

allocation_long <- allocations |>
  dplyr::select(
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
    y = value
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
  facet_wrap(
    ~ allocation_measure,
    scales = "free_y",
    ncol = 1
  ) +
  labs(
    title = "Joint allocation stability across random seeds",
    subtitle = "Narrow boxes indicate similar allocations across runs",
    x = "Stratum",
    y = "Allocated value"
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
# 8. Feasibility and repair diagnostics across seeds
# ============================================================

repair_diagnostics <- run_summary |>
  dplyr::select(
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
    y = value
  )
) +
  geom_line(aes(group = diagnostic), alpha = 0.6) +
  geom_point(size = 2) +
  facet_wrap(
    ~ diagnostic,
    scales = "free_y",
    ncol = 1
  ) +
  labs(
    title = "Repair behaviour across random seeds",
    x = "Random seed",
    y = "Diagnostic value"
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
