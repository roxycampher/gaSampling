repair_ckh <- function(c_h,
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
                       cost_unit = NULL) {
  L <- length(c_h)
  use_cost <- !is.null(max_cost)
  
  if (is.null(lower_c)) lower_c <- rep(2L, L)
  if (is.null(upper_c)) upper_c <- rep(Inf, L)
  if (is.null(lower_k)) lower_k <- rep(1L, L)
  if (is.null(upper_k)) upper_k <- rep(Inf, L)
  
  if (use_cost) {
    if (is.null(cost_cluster) || is.null(cost_unit)) {
      stop("cost_cluster and cost_unit must be supplied when max_cost is used.")
    }
    if (length(cost_cluster) == 1L) cost_cluster <- rep(cost_cluster, L)
    if (length(cost_unit) == 1L) cost_unit <- rep(cost_unit, L)
  }
  
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
    
    list(
      c_h = as.integer(c_candidate),
      k_h = as.integer(k_candidate),
      n_h = n_h,
      total_n = sum(n_h),
      target_total_n = total_n,
      total_cost = total_cost,
      max_cost = max_cost,
      cost_gap = if (use_cost) total_cost - max_cost else NA_real_,
      cost_ok = if (use_cost) total_cost <= max_cost else TRUE,
      sample_size_gap = sum(n_h) - total_n,
      repaired_exactly = sum(n_h) == total_n,
      iterations = iterations
    )
  }
  
  score_result <- function(result) {
    c_change <- sum(abs(result$c_h - start_c_h))
    k_change <- sum(abs(result$k_h - start_k_h))
    cost_over <- if (use_cost) max(0, result$cost_gap) else 0
    
    c(
      abs(result$sample_size_gap),
      cost_over,
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
      
      if (diff == 0L) {
        best_result <- current_result
        break
      }
      
      if (diff > 0) {
        possible <- which(k_candidate < upper_k & c_candidate <= diff)
        
        if (length(possible) > 0L) {
          j <- possible[which.max(c_candidate[possible])]
          jump <- min(
            floor(diff / c_candidate[j]),
            upper_k[j] - k_candidate[j],
            max(1L, as.integer(k_tolerance))
          )
          k_candidate[j] <- k_candidate[j] + as.integer(jump)
        } else {
          possible <- which(k_candidate < upper_k)
          if (length(possible) == 0L) break
          
          j <- possible[which.min(c_candidate[possible])]
          trial_k <- k_candidate
          trial_k[j] <- trial_k[j] + 1L
          
          if (abs(sum(c_candidate * trial_k) - total_n) >=
              abs(best_result$sample_size_gap)) {
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

repair_joint <- function(...) {
  repair_ckh(...)
}
