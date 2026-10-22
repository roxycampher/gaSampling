#' Joint two-stage allocation repair function
#'
#' Adjust cluster and unit allocations to meet sampling constraints
#'
#' Search for an integer allocation that meets or exceeds a minimum
#' total sample size, respects stratum allocation bounds and satisfies
#' an optional budget range.
#'
#' @param c_h Numeric vector giving the initial number of primary
#'   sampling units (clusters) sampled in each stratum.
#' @param k_h Numeric vector giving the initial number of secondary
#'   sampling units sampled per cluster in each stratum.
#' @param total_n A single non-negative integer giving the minimum
#'   required total sample size.
#' @param lower_c Numeric vector of minimum cluster counts per stratum.
#'   Defaults to two per stratum when NULL.
#' @param upper_c Numeric vector of maximum cluster counts per stratum.
#'   Defaults to Inf per stratum when NULL.
#' @param lower_k Numeric vector of minimum numbers of sampled units
#'   per cluster in each stratum. Defaults to one per stratum when NULL.
#' @param upper_k Numeric vector of maximum numbers of sampled units
#'   per cluster in each stratum. Defaults to Inf per stratum when NULL.
#' @param max_iter Non-negative integer giving the maximum number of
#'   unit-allocation adjustment steps for each candidate cluster
#'   allocation. Defaults to 1000.
#' @param k_tolerance Positive integer limiting the change in k_h
#'   within a single adjustment step. Defaults to 5.
#' @param c_search_radius Non-negative integer giving the maximum
#'   change explored in an individual stratum's initial cluster count.
#'   Also limits the number of successive cluster reductions considered
#'   when the initial allocation exceeds max_cost. Defaults to 2.
#' @param max_cost Positive numeric value giving the reference budget.
#'   If NULL, budget constraints are not applied.
#' @param cost_cluster Numeric vector of costs per sampled cluster in
#'   each stratum, or a single value applied to all strata. Required
#'   when max_cost is supplied.
#' @param cost_unit Numeric vector of costs per sampled secondary unit
#'   in each stratum, or a single value applied to all strata. Values
#'   should be positive. Required when max_cost is supplied.
#' @param min_budget_use Numeric proportion giving the minimum required
#'   budget use. Defaults to 0.95, requiring at least 95 percent of
#'   max_cost to be used. NULL is treated as zero.
#' @param budget_tol Non-negative numeric proportion allowing spending
#'   above max_cost. Defaults to 0.02, allowing spending up to 102
#'   percent of max_cost. NULL is treated as zero.
#'
#' @return A named list containing the selected allocation and diagnostics:
#'   \item{c_h}{Integer vector of selected cluster counts by stratum.}
#'   \item{k_h}{Integer vector of selected unit counts per cluster
#'     by stratum.}
#'   \item{n_h}{Integer vector of stratum sample sizes, c_h * k_h.}
#'   \item{total_n}{Total sample size achieved by the selected allocation.}
#'   \item{target_total_n}{The supplied minimum sample size requirement.}
#'   \item{total_cost}{Cost of the selected allocation, or NA when
#'     budget constraints are not applied.}
#'   \item{max_cost}{The supplied reference budget, or NULL.}
#'   \item{lower_cost_bound, upper_cost_bound}{Allowed spending limits,
#'     or NA when budget constraints are not applied.}
#'   \item{cost_gap}{Total cost minus max_cost, or NA when budget
#'     constraints are not applied.}
#'   \item{cost_over}{Amount above the upper spending limit, or zero
#'     when that limit is met or budget constraints are not applied.}
#'   \item{budget_underuse}{Amount below the lower spending limit, or
#'     zero when that limit is met or budget constraints are not applied.}
#'   \item{budget_use_pct}{Percentage of max_cost used, or NA when
#'     budget constraints are not applied.}
#'   \item{lower_cost_ok, upper_cost_ok, cost_ok}{Logical indicators
#'     for meeting the lower, upper, and combined budget constraints.
#'     All are TRUE when budget constraints are not applied.}
#'   \item{sample_size_gap}{Achieved sample size minus target_total_n.}
#'   \item{sample_shortfall}{Amount by which the achieved sample size
#'     falls below target_total_n, or zero if the requirement is met.}
#'   \item{sample_size_ok}{Logical value indicating whether the minimum
#'     sample size is met or exceeded.}
#'   \item{repaired_exactly}{Alias for sample_size_ok. Does not indicate
#'     exact equality with the target or satisfaction of the budget.}
#'   \item{iterations}{Iteration index recorded for the selected result
#'     within its unit-allocation search, rather than the total number
#'     of iterations across all candidates.}
#'
#' @details
#'   Allocation and bound vectors must have the same length, with
#'   entries ordered consistently by stratum. Scalar cost inputs are
#'   repeated across strata. Bounds should be positive integer-valued
#'   counts, with each lower bound no greater than its upper bound.
#'   Upper bounds may be Inf.
#'
#'   Initial allocations are rounded and restricted to their bounds.
#'   The function then adjusts unit counts for the initial cluster
#'   allocation and nearby candidate cluster allocations.
#'
#'   Candidates are compared in priority order by sample shortfall,
#'   spending above the upper budget limit, spending below the lower
#'   budget limit, absolute distance from max_cost and changes in
#'   cluster and unit counts from the rounded, bounded starting values.
#'   The search stops when an allocation satisfies both the sample
#'   size and budget requirements.
#'
#'   This is a heuristic search and does not guarantee that a feasible
#'   allocation will be found even when one exists. Check sample_size_ok
#'   and cost_ok before using the returned allocation.
#'   The function requires cost_twostage when max_cost is supplied.
#'   It checks vector lengths and required cost inputs, but does not
#'   fully validate input values or bound consistency.
#'
#' @examples
#' repair_joint(
#'   c_h = c(10, 15),
#'   k_h = c(4, 5),
#'   total_n = 150,
#'   lower_c = c(2, 2),
#'   upper_c = c(20, 25),
#'   lower_k = c(1, 1),
#'   upper_k = c(10, 10)
#' )
#' @export
repair_joint <- function(c_h,
                         k_h,
                         total_n,
                         lower_c = NULL,
                         upper_c = NULL,
                         lower_k = NULL,
                         upper_k = NULL,
                         max_iter = 1000,
                         k_tolerance = 5L,
                         c_search_radius = 2L,
                         max_cost = NULL,
                         cost_cluster = NULL,
                         cost_unit = NULL,
                         min_budget_use = 0.95,
                         budget_tol = 0.02) {
  L <- length(c_h)
  use_cost <- !is.null(max_cost)
  
  # Default
  if (is.null(lower_c)) lower_c <- rep(2L, L)
  if (is.null(upper_c)) upper_c <- rep(Inf, L)
  if (is.null(lower_k)) lower_k <- rep(1L, L)
  if (is.null(upper_k)) upper_k <- rep(Inf, L)
  
  if (use_cost) {
    if (is.null(cost_cluster) || is.null(cost_unit)) {
      stop("cost_cluster and cost_unit must be supplied when max_cost is used.")
    }
    if (is.null(min_budget_use)) min_budget_use <- 0
    if (is.null(budget_tol)) budget_tol <- 0
    if (length(cost_cluster) == 1L) cost_cluster <- rep(cost_cluster, L)
    if (length(cost_unit) == 1L) cost_unit <- rep(cost_unit, L)
  }
  
  # Length of vectors check
  if (length(k_h) != L ||
      length(lower_c) != L ||
      length(upper_c) != L ||
      length(lower_k) != L ||
      length(upper_k) != L ||
      (use_cost && length(cost_cluster) != L) ||
      (use_cost && length(cost_unit) != L)) {
    stop("c_h, k_h, bounds, and costs must have the same length.")
  }
  
  c_h <- as.integer(pmax(lower_c, pmin(upper_c, round(c_h))))
  k_h <- as.integer(pmax(lower_k, pmin(upper_k, round(k_h))))
  
  start_c_h <- c_h
  start_k_h <- k_h
  
  make_result <- function(c_candidate, k_candidate, iterations) {
    n_h <- as.integer(c_candidate * k_candidate)
    achieved_total_n <- sum(n_h)
    total_cost <- if (use_cost) {
      cost_twostage(
        c_h = c_candidate,
        k_h = k_candidate,
        cost_cluster = cost_cluster,
        cost_unit = cost_unit
      )
    } else {
      NA_real_
    }
    sample_size_gap <- achieved_total_n - total_n
    sample_shortfall <- max(0, total_n - achieved_total_n)
    sample_size_ok <- achieved_total_n >= total_n
    
    lower_cost_bound <- if (use_cost) max_cost * min_budget_use else NA_real_
    upper_cost_bound <- if (use_cost) max_cost * (1 + budget_tol) else NA_real_
    lower_cost_ok <- if (use_cost) total_cost >= lower_cost_bound else TRUE
    upper_cost_ok <- if (use_cost) total_cost <= upper_cost_bound else TRUE
    cost_ok <- lower_cost_ok && upper_cost_ok
    budget_underuse <- if (use_cost) max(0, lower_cost_bound - total_cost) else 0
    cost_over <- if (use_cost) max(0, total_cost - upper_cost_bound) else 0
    
    list(
      c_h = as.integer(c_candidate),
      k_h = as.integer(k_candidate),
      n_h = n_h,
      total_n = achieved_total_n,
      target_total_n = total_n,
      total_cost = total_cost,
      max_cost = max_cost,
      lower_cost_bound = lower_cost_bound,
      upper_cost_bound = upper_cost_bound,
      cost_gap = if (use_cost) total_cost - max_cost else NA_real_,
      cost_over = cost_over,
      budget_underuse = budget_underuse,
      budget_use_pct = if (use_cost) 100 * total_cost / max_cost else NA_real_,
      lower_cost_ok = lower_cost_ok,
      upper_cost_ok = upper_cost_ok,
      cost_ok = cost_ok,
      sample_size_gap = achieved_total_n - total_n,
      sample_shortfall = sample_shortfall,
      sample_size_ok = sample_size_ok,
      repaired_exactly = sample_size_ok,
      iterations = iterations
    )
  }
  
  score_result <- function(result) {
    c_change <- sum(abs(result$c_h - start_c_h))
    k_change <- sum(abs(result$k_h - start_k_h))
    unused_budget <- if (use_cost) abs(max_cost - result$total_cost) else 0
    
    c(
      result$sample_shortfall,
      result$cost_over,
      result$budget_underuse,
      unused_budget,
      c_change,
      k_change
    )
  }
  
  better_result <- function(candidate, current) {
    if (is.null(current)) return(TRUE)
    
    candidate_score <- score_result(candidate)
    current_score <- score_result(current)
    
    for (i in seq_along(candidate_score)) {
      if (candidate_score[i] < current_score[i]) return(TRUE)
      if (candidate_score[i] > current_score[i]) return(FALSE)
    }
    
    FALSE
  }
  
  repair_k_for_c <- function(c_candidate, k_reference) {
    c_candidate <- as.integer(pmax(lower_c, pmin(upper_c, c_candidate)))
    k_candidate <- as.integer(pmax(lower_k, pmin(upper_k, k_reference)))
    best_result <- make_result(c_candidate, k_candidate, 0L)
    
    for (step in seq_len(max_iter)) {
      current_total <- sum(c_candidate * k_candidate)
      diff <- total_n - current_total
      current_result <- make_result(c_candidate, k_candidate, step)
      
      if (better_result(current_result, best_result)) {
        best_result <- current_result
      }
      
      if (current_result$sample_size_ok && current_result$cost_ok) {
        best_result <- current_result
        break
      }
      
      if (diff > 0 || (use_cost && current_result$budget_underuse > 0)) {
        budget_room <- if (use_cost) {
          floor((max_cost * (1 + budget_tol) - current_result$total_cost) /
                  (cost_unit * c_candidate))
        } else {
          rep(Inf, length(k_candidate))
        }
        possible <- which(
          k_candidate < upper_k &
            budget_room > 0 &
            (diff > 0 | current_result$budget_underuse > 0)
        )
        
        if (diff > 0) {
          possible <- possible[c_candidate[possible] <= diff]
        }
        
        if (length(possible) == 0L && use_cost && current_result$budget_underuse > 0) {
          possible <- which(k_candidate < upper_k & budget_room > 0)
        }
        
        if (length(possible) > 0L) {
          j <- possible[which.max(c_candidate[possible])]
          jump <- min(
            if (diff > 0) max(1L, floor(diff / c_candidate[j])) else Inf,
            upper_k[j] - k_candidate[j],
            budget_room[j],
            max(1L, as.integer(k_tolerance))
          )
          k_candidate[j] <- k_candidate[j] + as.integer(jump)
        } else {
          possible <- which(k_candidate < upper_k)
          if (length(possible) == 0L) break
          
          j <- possible[which.min(c_candidate[possible])]
          trial_k <- k_candidate
          trial_k[j] <- trial_k[j] + 1L
          
          if (max(0, total_n - sum(c_candidate * trial_k)) >=
              best_result$sample_shortfall) {
            break
          }
          
          k_candidate <- trial_k
        }
      } else {
        diff_abs <- abs(diff)
        possible <- which(k_candidate > lower_k & c_candidate <= diff_abs)
        
        if (length(possible) > 0L) {
          j <- possible[which.max(c_candidate[possible])]
          jump <- min(
            floor(diff_abs / c_candidate[j]),
            k_candidate[j] - lower_k[j],
            max(1L, as.integer(k_tolerance))
          )
          k_candidate[j] <- k_candidate[j] - as.integer(jump)
        } else {
          possible <- which(k_candidate > lower_k)
          if (length(possible) == 0L) break
          
          j <- possible[which.min(c_candidate[possible])]
          trial_k <- k_candidate
          trial_k[j] <- trial_k[j] - 1L
          
          if (abs(sum(c_candidate * trial_k) - total_n) >=
              abs(best_result$sample_size_gap)) {
            break
          }
          
          k_candidate <- trial_k
        }
      }
    }
    
    best_result
  }
  
  c_candidates <- list(start_c_h)
  
  for (j in seq_len(L)) {
    for (delta in seq.int(-c_search_radius, c_search_radius)) {
      if (delta == 0L) next
      
      candidate <- start_c_h
      candidate[j] <- candidate[j] + delta
      candidate <- as.integer(pmax(lower_c, pmin(upper_c, candidate)))
      c_candidates[[length(c_candidates) + 1L]] <- candidate
    }
  }
  
  if (use_cost) {
    candidate_cost <- cost_twostage(
      c_h = start_c_h,
      k_h = start_k_h,
      cost_cluster = cost_cluster,
      cost_unit = cost_unit
    )
    
    if (candidate_cost > max_cost) {
      reduce_candidate <- start_c_h
      reductions <- 0L
      
      while (reductions < c_search_radius && sum(reduce_candidate) > sum(lower_c)) {
        reducible <- which(reduce_candidate > lower_c)
        if (length(reducible) == 0L) break
        
        saving <- cost_cluster[reducible] +
          cost_unit[reducible] * start_k_h[reducible]
        j <- reducible[which.max(saving)]
        reduce_candidate[j] <- reduce_candidate[j] - 1L
        reductions <- reductions + 1L
        
        c_candidates[[length(c_candidates) + 1L]] <- reduce_candidate
      }
    }
  }
  
  c_keys <- vapply(c_candidates, paste, collapse = ":", character(1))
  c_candidates <- c_candidates[!duplicated(c_keys)]
  
  best <- NULL
  
  for (candidate in c_candidates) {
    candidate_result <- repair_k_for_c(candidate, start_k_h)
    
    if (better_result(candidate_result, best)) {
      best <- candidate_result
    }
    
    if (best$repaired_exactly && best$cost_ok) {
      break
    }
  }
  
  best
}