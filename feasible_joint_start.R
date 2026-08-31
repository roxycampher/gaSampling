# =============================================================================
# gaSampling: feasible starting allocation for the joint c_h, k_h GA
#
# Finds a total number of sampled clusters and repaired starting allocation
# that satisfy the joint GA feasibility constraints.
# =============================================================================

allocate_ch_by_weight <- function(weights,
                                  total_c,
                                  upper,
                                  lower = rep(2L, length(weights))) {
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

find_feasible_joint_start <- function(strata,
                                      total_n,
                                      max_cost,
                                      max_candidate_c = 1000L,
                                      candidate_c = NULL,
                                      lower_c = NULL,
                                      upper_c = NULL,
                                      lower_k = NULL,
                                      upper_k = NULL,
                                      default_starting_suggestion =
                                        "3. Neyman proportional",
                                      require_all_starts = TRUE) {
  required_cols <- c(
    "stratum",
    "C_h",
    "N_h",
    "sigma2_bh",
    "cost_cluster",
    "cost_unit"
  )
  missing_cols <- setdiff(required_cols, names(strata))

  if (length(missing_cols) > 0L) {
    stop(
      "strata is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  L <- nrow(strata)
  
  # defaults
  if (is.null(lower_c)) lower_c <- rep(2L, L)
  if (is.null(upper_c)) upper_c <- strata$C_h
  if (is.null(lower_k)) lower_k <- rep(1L, L)
  if (is.null(upper_k)) upper_k <- rep(total_n, L)
  
  # The different suggestions
  starting_cases <- list(
    "1. Proportional to stratum size" = strata$N_h,
    "2. Proportional to between-cluster SD" = sqrt(strata$sigma2_bh),
    "3. Neyman proportional" = strata$N_h * sqrt(strata$sigma2_bh),
    "4. Equal c_h allocation" = rep(2, L)
  )

  if (!default_starting_suggestion %in% names(starting_cases)) {
    stop("default_starting_suggestion must be one of: ",
         paste(names(starting_cases), collapse = ", "))
  }

  if (is.null(candidate_c)) {
    lower_total_c <- sum(lower_c)
    upper_total_c <- min(
      sum(upper_c) - 2L,
      total_n,
      as.integer(max_candidate_c)
    )

    if (upper_total_c < lower_total_c) {
      stop("No non-saturated total_c is available for the starting search.")
    }

    candidate_c <- unique(as.integer(c(
      round(seq(upper_total_c, lower_total_c, length.out = 60)),
      90L,
      250L,
      500L
    )))

    candidate_c <- sort(
      candidate_c[
        candidate_c >= lower_total_c &
          candidate_c <= upper_total_c
      ],
      decreasing = TRUE
    )
  } else {
    candidate_c <- sort(unique(as.integer(candidate_c)), decreasing = TRUE)
  }

  attempted <- list()

  for (current_total_c in candidate_c) {
    starts <- list()

    for (starting_name in names(starting_cases)) {
      initial_c_h <- allocate_ch_by_weight(
        weights = starting_cases[[starting_name]],
        total_c = current_total_c,
        upper = upper_c,
        lower = lower_c
      )

      initial_k_h <- rep( # set the k_h to equal for the c_h allocation suggestions
        max(1L, floor(total_n / sum(initial_c_h))),
        L
      )

      repaired <- repair_ckh(
        c_h = initial_c_h,
        k_h = initial_k_h,
        total_n = total_n,
        lower_c = lower_c,
        upper_c = upper_c,
        lower_k = lower_k,
        upper_k = upper_k
      )

      constraint_check <- check_constraints(
        c_h = repaired$c_h,
        k_h = repaired$k_h,
        total_n = total_n,
        max_cost = max_cost,
        cost_cluster = strata$cost_cluster,
        cost_unit = strata$cost_unit,
        min_clusters = lower_c,
        min_units = lower_k,
        max_clusters = upper_c,
        max_units = upper_k
      )

      starts[[starting_name]] <- list(
        starting_suggestion = starting_name,
        initial_c_h = initial_c_h,
        initial_k_h = initial_k_h,
        repair = repaired,
        constraint_check = constraint_check,
        feasible = constraint_check$feasible && repaired$repaired_exactly
      )
    }

    feasible_by_start <- vapply(
      starts,
      function(x) isTRUE(x$feasible),
      logical(1)
    )

    attempted[[as.character(current_total_c)]] <- feasible_by_start

    if (
      (require_all_starts && all(feasible_by_start)) ||
        (!require_all_starts && feasible_by_start[[default_starting_suggestion]])
    ) {
      selected <- starts[[default_starting_suggestion]]

      all_starts <- do.call(
        rbind,
        lapply(starts, function(x) {
          data.frame(
            starting_suggestion = x$starting_suggestion,
            stratum = strata$stratum,
            initial_c_h = x$initial_c_h,
            initial_k_h = x$initial_k_h,
            repaired_c_h = x$repair$c_h,
            repaired_k_h = x$repair$k_h,
            repaired_n_h = x$repair$c_h * x$repair$k_h,
            repaired_total_n = x$repair$total_n,
            repaired_total_c = sum(x$repair$c_h),
            repaired_total_cost = x$constraint_check$total_cost,
            sample_size_gap = x$repair$sample_size_gap,
            repaired_exactly = x$repair$repaired_exactly,
            feasible = x$feasible,
            repair_iterations = x$repair$iterations,
            row.names = NULL
          )
        })
      )

      diagnostics <- aggregate(
        cbind(
          repaired_total_n,
          repaired_total_c,
          repaired_total_cost,
          sample_size_gap,
          repair_iterations
        ) ~ starting_suggestion,
        data = all_starts,
        FUN = function(x) x[1]
      )
      diagnostics$repaired_exactly <- aggregate(
        repaired_exactly ~ starting_suggestion,
        data = all_starts,
        FUN = function(x) x[1]
      )$repaired_exactly
      diagnostics$feasible <- aggregate(
        feasible ~ starting_suggestion,
        data = all_starts,
        FUN = function(x) x[1]
      )$feasible

      return(list(
        total_c = current_total_c,
        c_h = selected$repair$c_h,
        k_h = selected$repair$k_h,
        suggestion = c(selected$repair$c_h, selected$repair$k_h),
        starting_suggestion = default_starting_suggestion,
        all_starts = all_starts,
        diagnostics = diagnostics,
        attempted = attempted
      ))
    }
  }

  stop("Could not find a feasible joint starting allocation.")
}
