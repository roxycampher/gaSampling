rm(list = ls())

library(dplyr)
library(readr)
library(R2BEAT)

source("variance_twostage.R")
source("cost_twostage.R")
source("read_survey_data.R")
source("sampling_frame.R")
source("repair_joint.R")

################################################################################
# 1. READ DATA
################################################################################

toy_data <- read.csv("ecdc-2021-public-v2.csv")

toy_data$count_children_present <- readr::parse_number(
  toy_data$count_children_present
)

# Remove observations without a usable target value
r2b_frame <- toy_data |>
  filter(
    !is.na(province),
    !is.na(id_ward),
    !is.na(count_children_present)
  )

################################################################################
# 2. CREATE R2BEAT-SAFE IDENTIFIERS
################################################################################

# Numeric stratum identifier
province_lookup <- r2b_frame |>
  distinct(province) |>
  arrange(province) |>
  mutate(R2B_STRATUM = row_number())

r2b_frame <- r2b_frame |>
  left_join(
    province_lookup,
    by = "province"
  )

# PSU must be unique across the whole frame.
# interaction(province, id_ward) protects us in case ward IDs are only
# unique within province

r2b_frame <- r2b_frame |>
  mutate(
    PSU_KEY = interaction(
      province,
      id_ward,
      drop = TRUE
    ),
    R2B_PSU = as.integer(PSU_KEY),
    
    # Unique final sampling unit ID
    R2B_SSU = row_number(),
    
    # One national domain
    R2B_DOMAIN = 1L,
    
    # Used internally by R2BEAT
    one = 1
  )

################################################################################
# 3. CHECK THE FRAME
################################################################################

cat("\nNumber of observations:", nrow(r2b_frame), "\n")
cat("Number of strata:", n_distinct(r2b_frame$R2B_STRATUM), "\n")
cat("Number of PSUs:", n_distinct(r2b_frame$R2B_PSU), "\n")

################################################################################
# 4. PREPARE R2BEAT INPUTS
################################################################################

# minimum:
# Minimum SSUs selected per selected PSU.
#
# Setting this to 1 makes the R2BEAT constraint as close as possible to
# k_h >= 1 in our GA
#
# Later run a sensitivity test with larger values such as 5 or 10 (check)

minimum_r2b <- 1

# Each row is one final survey unit
delta_r2b <- 1

# Suggested DEFF used when preparing the input.
#
# Note: when beat.2st() is later called with deft_start = NULL,
# the iterative allocation starts at DEFT = 1 and computes its own
# design-effect corrections

deff_sugg_r2b <- 1

r2b_input <- prepareInputToAllocation1(
  samp_frame  = r2b_frame,
  id_PSU      = "R2B_PSU",
  id_SSU      = "R2B_SSU",
  strata_var  = "R2B_STRATUM",
  target_vars = "count_children_present",
  deff_var    = "R2B_STRATUM",
  domain_var  = "R2B_DOMAIN",
  minimum     = minimum_r2b,
  delta       = delta_r2b,
  deff_sugg   = deff_sugg_r2b
)

################################################################################
# 5. INSPECT THE GENERATED INPUTS
################################################################################

cat("\n--- R2BEAT STRATA INPUT ---\n")
print(r2b_input$strata)

cat("\n--- R2BEAT ICC / RHO ---\n")
print(r2b_input$rho)

cat("\n--- R2BEAT DESIGN FILE ---\n")
print(r2b_input$des_file)

cat("\n--- R2BEAT PSU FILE ---\n")
print(head(r2b_input$psu_file))

################################################################################
# 6. CHECK INPUT COHERENCE
################################################################################

check_input(
  strata = r2b_input$strata,
  des = r2b_input$des_file,
  strata_var_strata = "STRATUM",
  strata_var_des = "STRATUM"
)

################################################################################
# 7. FUNCTION TO RUN R2BEAT FOR A GIVEN CV
################################################################################

run_r2beat <- function(cv_target,
                       input,
                       min_PSU_strat = 2) {
  
  # One target variable, one national domain
  errors <- data.frame(
    DOM = "DOM1",
    CV1 = cv_target
  )
  
  result <- beat.2st(
    stratif = input$strata,
    errors = errors,
    des_file = input$des_file,
    psu_file = input$psu_file,
    rho = input$rho,
    
    # Start at DEFT = 1.
    # R2BEAT subsequently estimates design effects iteratively.
    deft_start = NULL,
    
    effst = input$effst,
    
    minPSUstrat = min_PSU_strat,
    minnumstrat = 2
  )
  
  result
}

################################################################################
# 8. CALIBRATE CV SO R2BEAT RETURNS APPROXIMATELY total_n
################################################################################

target_n <- 100000

# Smaller CV -> greater precision -> larger sample
# Larger CV  -> lower precision   -> smaller sample

cv_grid <- seq(
  from = 0.001,
  to = 0.02,
  by = 0.0005
)

calibration_results <- data.frame(
  CV = cv_grid,
  total_n = NA_real_
)

r2beat_runs <- vector(
  "list",
  length(cv_grid)
)

for (i in seq_along(cv_grid)) {
  
  cat(
    "\nRunning CV =",
    cv_grid[i],
    "\n"
  )
  
  fit <- try(
    run_r2beat(
      cv_target = cv_grid[i],
      input = r2b_input
    ),
    silent = TRUE
  )
  
  if (!inherits(fit, "try-error")) {
    
    r2beat_runs[[i]] <- fit
    
    # Final SSU allocation
    allocation_table <- fit$alloc
    
    # Remove the "Total" row
    allocation_table <- allocation_table[
      allocation_table$STRATUM != "Total",
    ]
    
    calibration_results$total_n[i] <- sum(
      allocation_table$ALLOC
    )
  }
}

print(calibration_results)

################################################################################
# 9. FIND CV GIVING SAMPLE SIZE CLOSEST TO N_TOTAL
################################################################################

calibration_results <- calibration_results |>
  mutate(
    difference_from_target = abs(total_n - target_n)
  )

best_index <- which.min(
  calibration_results$difference_from_target
)

best_cv <- calibration_results$CV[best_index]
#best_cv <- 0.001

cat(
  "\nBest CV:",
  best_cv,
  "\n"
)

cat(
  "R2BEAT sample size:",
  calibration_results$total_n[best_index],
  "\n"
)

r2beat_result <- r2beat_runs[[best_index]]

print(r2beat_result$iterations)
print(r2beat_result$alloc)
print(r2beat_result$expected)

################################################################################
# 10. EXTRACT STRATUM SSU ALLOCATION n_h
################################################################################

r2b_nh <- r2beat_result$alloc |>
  filter(
    STRATUM != "Total"
  ) |>
  transmute(
    R2B_STRATUM = as.integer(
      as.character(STRATUM)
    ),
    n_h_r2beat = as.integer(ALLOC)
  )

r2b_nh <- r2b_nh |>
  left_join(
    province_lookup,
    by = "R2B_STRATUM"
  )

print(r2b_nh)

################################################################################
# 11. EXTRACT NUMBER OF PSUs c_h
################################################################################

psu_structure <- r2beat_result$psu_trs

min_PSU_strat <- 2

r2b_ch <- lapply(
  split(
    psu_structure,
    psu_structure$STRATUM
  ),
  function(x) {
    
    # Self-representative PSUs
    n_SR <- sum(
      x$AR == 1
    )
    
    # Number of non-self-representative PSU substrata
    nsr_substrata <- unique(
      x$SUBSTRAT[x$AR == 0]
    )
    
    n_NSR <- length(
      nsr_substrata
    ) * min_PSU_strat
    
    data.frame(
      R2B_STRATUM = as.integer(
        as.character(x$STRATUM[1])
      ),
      c_h_r2beat = n_SR + n_NSR,
      PSU_SR = n_SR,
      PSU_NSR = n_NSR
    )
  }
) |>
  bind_rows()

r2b_ch <- r2b_ch |>
  left_join(
    province_lookup,
    by = "R2B_STRATUM"
  )

print(r2b_ch)

################################################################################
# 12. CREATE R2BEAT c_h / k_h ALLOCATION
################################################################################

r2beat_allocation_original <- r2b_nh |>
  select(
    R2B_STRATUM,
    province,
    n_h_r2beat
  ) |>
  left_join(
    r2b_ch |>
      select(
        R2B_STRATUM,
        c_h_r2beat,
        PSU_SR,
        PSU_NSR
      ),
    by = "R2B_STRATUM"
  ) |>
  mutate(
    
    # Effective average SSUs sampled per PSU
    k_h_r2beat = n_h_r2beat / c_h_r2beat
    
  )

print(r2beat_allocation_original)

# we need to compare the k_h integer solution as well - but isn't it bad to make it cost aware?
L <- nrow(r2beat_allocation_original)

lower_c <- rep(1L, L)
upper_c <- rep(n_distinct(r2b_frame$R2B_PSU), L)

lower_k <- rep(1L, L)
upper_k <- rep(target_n, L)

r2beat_allocation_cost_aware <- repair_ckh(c_h = r2beat_allocation_original$c_h_r2beat,
                                         k_h = r2beat_allocation_original$k_h_r2beat,
                                         total_n = target_n,
                                         lower_c = lower_c,
                                         upper_c = upper_c,
                                         lower_k = lower_k,
                                         upper_k = upper_k,
                                         max_iter = 1000,
                                         k_tolerance = 5L,
                                         c_search_radius = 2L,
                                         max_cost = 1000000,
                                         cost_cluster = 1000,
                                         cost_unit = 10
                                         )

print(r2beat_allocation_cost_aware)

################################################################################
# 13. PREPARE STRATA DATA
################################################################################

survey_data <- read_survey_data(
  data = toy_data,
  stratum_id = "province",
  cluster_id = "id_ward",
  y_var = "count_children_present",
  na_action = "warn"
)

max_cost <- 1000000000
available_c <- nrow(unique(survey_data[c("strata", "cluster")]))
affordable_c <- floor((max_cost - 10 * target_n) / 100)
full_frame_c <- min(available_c, affordable_c)

full_frame <- prepare_sampling_frame(
  data = survey_data,
  y_var = "count_children_present",
  total_c = full_frame_c,
  cost_cluster = 1000,
  cost_unit = 100
)

strata <- full_frame$strata

strata

L <- nrow(strata)

lower_c <- rep(1L, L)
upper_c <- strata$C_h

lower_k <- rep(1L, L)
upper_k <- rep(target_n, L)

################################################################################
# 14. ALIGN R2BEAT WITH STRATA ORDER
################################################################################

r2beat_allocation_original <- r2beat_allocation_original |>
  arrange(
    match(
      province,
      strata$province
    )
  )

c_r2beat_og <- r2beat_allocation_original$c_h_r2beat

k_r2beat_og <- r2beat_allocation_original$k_h_r2beat

c_r2beat_costaw <- r2beat_allocation_cost_aware$c_h

k_r2beat_costaw <- r2beat_allocation_cost_aware$k_h

################################################################################
# 15. EVALUATE WITH OUR VARIANCE FUNCTION
################################################################################

variance_r2beat_noninteger <- variance_twostage(
  c_h = r2beat_allocation_original$c_h_r2beat,
  k_h = r2beat_allocation_original$k_h_r2beat,
  W_h = strata$W_h, # the weight of the stratum n_h/n_total (is the same even with R2BEAT)
  sigma2_bh = strata$sigma2_bh,
  sigma2_wh = strata$sigma2_wh
)

variance_r2beat_noninteger

variance_r2beat_cost_aware <- variance_twostage(
  c_h = c_r2beat_costaw,
  k_h = k_r2beat_costaw,
  W_h = strata$W_h, # the weight of the stratum n_h/n_total (is the same even with R2BEAT)
  sigma2_bh = strata$sigma2_bh,
  sigma2_wh = strata$sigma2_wh
)

variance_r2beat_cost_aware

################################################################################
# 16. OVERALL COMPARISON
################################################################################

cost_r2beat_og <- cost_twostage(c_r2beat_og,
                                k_r2beat_og,
                                cost_cluster = 1000,
                                cost_unit = 100)

cost_r2beat_ca <- cost_twostage(c_r2beat_costaw,
                                k_r2beat_costaw,
                                cost_cluster = 1000,
                                cost_unit = 100)

comparison <- data.frame(
  method = c(
    "R2BEAT non-integer",
    "R2BEAT integer cost aware"
  ),
  
  N_total =c(
    sum(r2beat_allocation_original$n_h_r2beat),
    sum(r2beat_allocation_cost_aware$n_h)
  ),
  
  total_PSU = c(
    sum(c_r2beat_og),
    sum(c_r2beat_costaw)
  ),
  
  total_SSU = c(
    sum(r2beat_allocation_original$n_h_r2beat),
    sum(r2beat_allocation_cost_aware$n_h)
  ),
  
  variance = c(
    variance_r2beat_noninteger,
    variance_r2beat_cost_aware
  ),
  
  cost = c(
    cost_r2beat_og,
    cost_r2beat_ca
  )
)

print(comparison)
