rm(list = ls())

library(dplyr)
library(ggplot2)
library(readr)
library(GA)
library(tidyr)
library(scales)

source("read_survey_data.R")
source("sampling_frame.R")
source("variance_twostage.R")
source("cost_twostage.R")
source("check_constraints.R")
source("repair_kh.R")
source("sample_size_feasibility.R")
source("fitness_kh.R")
source("allocate_twostage_ga_kh.R")

output_dir <- "kh_diagnostics_output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

allocate_ch_proportional <- function(C_h, total_c, lower = rep(1L, length(C_h))) {
  C_h <- as.integer(C_h)
  lower <- as.integer(lower)
  total_c <- as.integer(total_c)
  
  if (length(lower) != length(C_h)) {
    stop("lower must have the same length as C_h.")
  }
  
  if (total_c < sum(lower) || total_c > sum(C_h)) {
    stop("total_c is infeasible for these C_h bounds.")
  }
  
  c_h <- lower
  remaining <- total_c - sum(c_h)
  capacity <- C_h - c_h
  
  if (remaining > 0L) {
    raw_extra <- remaining * capacity / sum(capacity)
    extra <- pmin(floor(raw_extra), capacity)
    c_h <- c_h + extra
    
    while (sum(c_h) < total_c) {
      eligible <- which(c_h < C_h)
      frac <- raw_extra - floor(raw_extra)
      h <- eligible[order(frac[eligible], C_h[eligible] - c_h[eligible], decreasing = TRUE)][1]
      c_h[h] <- c_h[h] + 1L
    }
  }
  
  as.integer(c_h)
}

extract_ga_progress <- function(ga_result, case_name) {
  ga_summary <- as.data.frame(ga_result$ga_object@summary)
  ga_summary$generation <- seq_len(nrow(ga_summary))
  
  ga_summary |>
    transmute(
      case = case_name,
      generation = generation,
      best_fitness = max,
      mean_fitness = mean,
      q1_fitness = q1,
      median_fitness = median,
      q3_fitness = q3,
      min_fitness = min,
      best_variance = -best_fitness,
      mean_variance = -mean_fitness,
      q1_variance = -q3_fitness,
      median_variance = -median_fitness,
      q3_variance = -q1_fitness,
      worst_variance = -min_fitness,
      best_so_far_variance = cummin(best_variance)
    ) |>
    mutate(
      improvement_since_previous = lag(best_so_far_variance) - best_so_far_variance,
      new_best = best_so_far_variance < lag(best_so_far_variance, default = Inf)
    )
}

summarise_case <- function(case_name, allocation, strata, total_n, max_cost) {
  constraint_check <- check_constraints(
    c_h = allocation$c_h,
    k_h = allocation$k_h,
    total_n = total_n,
    max_cost = max_cost,
    cost_cluster = strata$cost_cluster,
    cost_unit = strata$cost_unit
  )
  
  data.frame(
    case = case_name,
    variance = variance_twostage(
      c_h = allocation$c_h,
      k_h = allocation$k_h,
      W_h = strata$W_h,
      sigma2_bh = strata$sigma2_bh,
      sigma2_wh = strata$sigma2_wh
    ),
    total_n = sum(allocation$n_h),
    total_c = sum(allocation$c_h),
    total_cost = cost_twostage(
      c_h = allocation$c_h,
      k_h = allocation$k_h,
      cost_cluster = strata$cost_cluster,
      cost_unit = strata$cost_unit
    ),
    max_cost = max_cost,
    feasible = constraint_check$feasible
  )
}

add_context <- function(allocation, case_name, strata) {
  allocation |>
    mutate(
      case = case_name,
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

toy_data <- read.csv("ecdc-2021-public-v2.csv")

clean_data <- read_survey_data(
  data = toy_data,
  stratum_id = "province",
  cluster_id = "id_ward",
  y_var = "count_children_present",
  na_action = "warn"
)

total_c <- 90
total_n <- 900
max_cost <- 20000

frame <- prepare_sampling_frame(
  data = clean_data,
  y_var = "count_children_present",
  total_c = total_c,
  cost_cluster = 100,
  cost_unit = 10
)

strata <- frame$strata

c_h_equal <- rep(total_c / nrow(strata), nrow(strata))

c_h_proportional <- allocate_ch_proportional(
  C_h = strata$C_h,
  total_c = total_c,
  lower = rep(1L, nrow(strata))
)

# 1. Equal c_h and optimised k_h
ga_equal_ch <- allocate_twostage_ga_kh(
  strata = strata,
  total_n = total_n,
  max_cost = max_cost,
  c_h = c_h_equal,
  popSize = 100,
  maxiter = 200,
  seed = 123,
  run = 200
)

case_1 <- ga_equal_ch$allocation

# 2. Proportional c_h and optimised k_h
ga_proportional_ch <- allocate_twostage_ga_kh(
  strata = strata,
  total_n = total_n,
  max_cost = max_cost,
  c_h = c_h_proportional,
  popSize = 100,
  maxiter = 200,
  seed = 123,
  run = 200
)

case_2 <- ga_proportional_ch$allocation

# 3. Equal allocation benchmark: equal c_h and equal k_h
k_h_equal <- rep(total_n / total_c, nrow(strata))

case_3 <- data.frame(
  stratum = strata$stratum,
  c_h = c_h_equal,
  k_h = k_h_equal,
  n_h = c_h_equal * k_h_equal
)

# 4. Proportional allocation benchmark: proportional c_h and proportional n_h
k_h_proportional <- repair_kh(
  x = (total_n * strata$W_h) / c_h_proportional,
  c_h = c_h_proportional,
  total_n = total_n,
  lower = rep(1L, nrow(strata)),
  upper = floor(total_n / c_h_proportional)
)

case_4 <- data.frame(
  stratum = strata$stratum,
  c_h = c_h_proportional,
  k_h = k_h_proportional,
  n_h = c_h_proportional * k_h_proportional
)

case_allocations <- bind_rows(
  add_context(case_1, "1. Equal c_h, optimised k_h", strata),
  add_context(case_2, "2. Proportional c_h, optimised k_h", strata),
  add_context(case_3, "3. Equal allocation benchmark", strata),
  add_context(case_4, "4. Proportional allocation benchmark", strata)
)

case_summary <- bind_rows(
  summarise_case("1. Equal c_h, optimised k_h", case_1, strata, total_n, max_cost),
  summarise_case("2. Proportional c_h, optimised k_h", case_2, strata, total_n, max_cost),
  summarise_case("3. Equal allocation benchmark", case_3, strata, total_n, max_cost),
  summarise_case("4. Proportional allocation benchmark", case_4, strata, total_n, max_cost)
) |>
  arrange(variance) |>
  mutate(
    rank = row_number(),
    improvement_vs_equal_benchmark_pct =
      100 *
      (
        variance[case == "3. Equal allocation benchmark"][1] -
          variance
      ) /
      variance[case == "3. Equal allocation benchmark"][1]
  )

generation_progress <- bind_rows(
  extract_ga_progress(ga_equal_ch, "1. Equal c_h, optimised k_h"),
  extract_ga_progress(ga_proportional_ch, "2. Proportional c_h, optimised k_h")
)

ch_allocation_comparison <- bind_rows(
  data.frame(
    case = "Equal c_h",
    stratum = strata$stratum,
    C_h = strata$C_h,
    W_h = strata$W_h,
    c_h = c_h_equal
  ),
  data.frame(
    case = "Proportional c_h",
    stratum = strata$stratum,
    C_h = strata$C_h,
    W_h = strata$W_h,
    c_h = c_h_proportional
  )
) |>
  mutate(
    cluster_population_share = C_h / sum(strata$C_h),
    sampled_cluster_share = c_h / sum(c_h),
    cluster_sampling_fraction = c_h / C_h
  )

readr::write_csv(case_summary, file.path(output_dir, "kh_four_case_summary.csv"))
readr::write_csv(case_allocations, file.path(output_dir, "kh_four_case_allocations.csv"))
readr::write_csv(generation_progress, file.path(output_dir, "kh_four_case_generation_progress.csv"))
readr::write_csv(ch_allocation_comparison, file.path(output_dir, "kh_ch_allocation_comparison.csv"))

p_case_comparison <- ggplot(case_summary, aes(x = reorder(case, variance), y = variance, fill = case)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  labs(
    title = "Four allocation cases",
    subtitle = "Lower variance is better",
    x = NULL,
    y = "Variance"
  ) +
  theme_minimal(base_size = 12)

p_generation <- ggplot(
  generation_progress,
  aes(x = generation, y = best_so_far_variance, colour = case)
) +
  geom_line(linewidth = 1) +
  labs(
    title = "GA convergence for optimised k_h cases",
    subtitle = "Case 1 fixes c_h equally; case 2 fixes c_h proportionally",
    x = "Generation",
    y = "Best-so-far variance",
    colour = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

p_population_spread <- ggplot(generation_progress, aes(x = generation)) +
  geom_ribbon(aes(ymin = q1_variance, ymax = q3_variance, fill = case), alpha = 0.25) +
  geom_line(aes(y = best_so_far_variance, colour = case), linewidth = 1) +
  facet_wrap(~ case) +
  labs(
    title = "Population spread by generation",
    subtitle = "Ribbon shows the interquartile range; line is best-so-far",
    x = "Generation",
    y = "Variance",
    fill = NULL,
    colour = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

p_final_allocation <- ggplot(case_allocations, aes(x = reorder(stratum, n_h), y = n_h, fill = case)) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(
    title = "Final n_h allocation by case",
    x = NULL,
    y = "n_h",
    fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

p_sample_share <- case_allocations |>
  select(case, stratum, n_share, population_share) |>
  pivot_longer(
    cols = c(n_share, population_share),
    names_to = "share_type",
    values_to = "share"
  ) |>
  mutate(
    share_type = recode(
      share_type,
      n_share = "Sample share",
      population_share = "Population share"
    )
  ) |>
  ggplot(aes(x = reorder(stratum, share), y = share, fill = share_type)) +
  geom_col(position = "dodge") +
  coord_flip() +
  facet_wrap(~ case) +
  scale_y_continuous(labels = percent) +
  labs(
    title = "Sample share versus population share",
    x = NULL,
    y = "Share",
    fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

p_ch <- ch_allocation_comparison |>
  select(case, stratum, cluster_population_share, sampled_cluster_share) |>
  pivot_longer(
    cols = c(cluster_population_share, sampled_cluster_share),
    names_to = "share_type",
    values_to = "share"
  ) |>
  mutate(
    share_type = recode(
      share_type,
      cluster_population_share = "Available cluster share",
      sampled_cluster_share = "Sampled cluster share"
    )
  ) |>
  ggplot(aes(x = reorder(stratum, share), y = share, fill = share_type)) +
  geom_col(position = "dodge") +
  coord_flip() +
  facet_wrap(~ case) +
  scale_y_continuous(labels = percent) +
  labs(
    title = "c_h allocation comparison",
    x = NULL,
    y = "Share of clusters",
    fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(output_dir, "kh_four_case_variance_comparison.png"), p_case_comparison, width = 9, height = 5, dpi = 300)
ggsave(file.path(output_dir, "kh_four_case_ga_convergence.png"), p_generation, width = 9, height = 5, dpi = 300)
ggsave(file.path(output_dir, "kh_four_case_ga_population_spread.png"), p_population_spread, width = 9, height = 5, dpi = 300)
ggsave(file.path(output_dir, "kh_four_case_final_allocation.png"), p_final_allocation, width = 10, height = 6, dpi = 300)
ggsave(file.path(output_dir, "kh_four_case_sample_share.png"), p_sample_share, width = 11, height = 7, dpi = 300)
ggsave(file.path(output_dir, "kh_ch_allocation_comparison.png"), p_ch, width = 10, height = 6, dpi = 300)

cat("\nFour-case summary:\n")
print(case_summary)

cat("\nFinal allocations:\n")
print(case_allocations |> select(case, stratum, c_h, k_h, n_h))

cat("\nSaved k_h diagnostics to:", normalizePath(output_dir), "\n")

################ MULTIPLE RUNS ########################################

seeds <- 1:20

run_ga_case <- function(c_h, case_name, seeds) {
  
  runs <- lapply(seeds, function(s) {
    
    result <- allocate_twostage_ga_kh(
      strata = strata,
      total_n = total_n,
      max_cost = max_cost,
      c_h = c_h,
      popSize = 100,
      maxiter = 200,
      run = 200,
      seed = s
    )
    
    allocation <- result$allocation
    
    constraint_check <- check_constraints(
      c_h = allocation$c_h,
      k_h = allocation$k_h,
      total_n = total_n,
      max_cost = max_cost,
      cost_cluster = strata$cost_cluster,
      cost_unit = strata$cost_unit
    )
    
    list(
      seed = s,
      result = result,
      allocation = allocation,
      variance = variance_twostage(
        c_h = allocation$c_h,
        k_h = allocation$k_h,
        W_h = strata$W_h,
        sigma2_bh = strata$sigma2_bh,
        sigma2_wh = strata$sigma2_wh
      ),
      total_cost = cost_twostage(
        c_h = allocation$c_h,
        k_h = allocation$k_h,
        cost_cluster = strata$cost_cluster,
        cost_unit = strata$cost_unit
      ),
      feasible = constraint_check$feasible
    )
  })
  
  run_summary <- bind_rows(lapply(runs, function(x) {
    data.frame(
      case = case_name,
      seed = x$seed,
      variance = x$variance,
      total_cost = x$total_cost,
      feasible = x$feasible
    )
  }))
  
  allocations <- bind_rows(lapply(runs, function(x) {
    x$allocation |>
      mutate(
        case = case_name,
        seed = x$seed
      )
  }))
  
  list(
    runs = runs,
    summary = run_summary,
    allocations = allocations
  )
}

multi_equal <- run_ga_case(
  c_h = c_h_equal,
  case_name = "Equal c_h, optimised k_h",
  seeds = seeds
)

multi_proportional <- run_ga_case(
  c_h = c_h_proportional,
  case_name = "Proportional c_h, optimised k_h",
  seeds = seeds
)


ga_run_summary <- bind_rows(
  multi_equal$summary,
  multi_proportional$summary
)

ga_allocations <- bind_rows(
  multi_equal$allocations,
  multi_proportional$allocations
)

ga_stability_summary <- ga_run_summary |>
  group_by(case) |>
  summarise(
    runs = n(),
    best_variance = min(variance),
    median_variance = median(variance),
    mean_variance = mean(variance),
    sd_variance = sd(variance),
    q1_variance = quantile(variance, 0.25),
    q3_variance = quantile(variance, 0.75),
    worst_variance = max(variance),
    feasibility_rate = mean(feasible),
    .groups = "drop"
  )

print(ga_stability_summary)

p_seed_variance <- ggplot(
  ga_run_summary,
  aes(
    x = seed,
    y = variance,
    colour = case
  )
) +
  geom_point(size = 3, alpha = 0.85) +
  scale_x_continuous(
    breaks = sort(unique(ga_run_summary$seed))
  ) +
  labs(
    title = "Final GA variance across random seeds",
    subtitle = "Point colour indicates the fixed cluster allocation",
    x = "Seed",
    y = "Final variance",
    colour = expression(c[h] ~ "allocation")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor.x = element_blank()
  )

ggsave(
  file.path(output_dir, "kh_variance_by_seed.png"),
  p_seed_variance,
  width = 9,
  height = 5,
  dpi = 300
)

# Do strata receive similar k_h across runs?
allocation_stability <- ga_allocations |>
  group_by(case, stratum) |>
  summarise(
    mean_k_h = mean(k_h),
    median_k_h = median(k_h),
    sd_k_h = sd(k_h),
    min_k_h = min(k_h),
    max_k_h = max(k_h),
    .groups = "drop"
  )


p_kh_stability <- ggplot(
  ga_allocations,
  aes(x = factor(stratum), y = k_h, fill = case)
) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  labs(
    title = "Stability of k_h allocations across GA runs",
    x = "Stratum",
    y = "k_h",
    fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(output_dir, "kh_ch_stability.png"), p_kh_stability, width = 10, height = 6, dpi = 300)
