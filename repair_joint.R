repair_ckh <- function(c_h,
                       k_h,
                       total_n,
                       lower_c = NULL,
                       upper_c = NULL,
                       lower_k = NULL,
                       upper_k = NULL,
                       max_iter = 1000,
                       k_tolerance = 5L,
                       c_search_radius = 2L) { # how far the function may search around each orgininal ch
  
  L <- length(c_h) # number of strata
  
  # Default bounds
  if (is.null(lower_c)) lower_c <- rep(2L, L) # Assumes lower bounds for ch are >= 2
  if (is.null(upper_c)) upper_c <- rep(Inf, L) # C_h
  if (is.null(lower_k)) lower_k <- rep(1L, L) # Assumes we have to sample at least one PSU
  if (is.null(upper_k)) upper_k <- rep(Inf, L)
  
  # check that all vectors have the same length
  if (length(k_h) != L ||
      length(lower_c) != L ||
      length(upper_c) != L ||
      length(lower_k) != L ||
      length(upper_k) != L) {
    stop("c_h, k_h, and all bounds must have the same length.")
  }
  
  # Round and constrain starting solution
  c_h <- as.integer(pmax(lower_c, pmin(upper_c, round(c_h))))
  k_h <- as.integer(pmax(lower_k, pmin(upper_k, round(k_h))))
  
  # initial solution (rounded)
  start_c_h <- c_h
  start_k_h <- k_h
  
  # Helper function
  # packages the allocation into a result list
  make_result <- function(c_candidate, k_candidate, checked) {
    n_h <- as.integer(c_candidate * k_candidate)
    
    list(
      c_h = as.integer(c_candidate),
      k_h = as.integer(k_candidate),
      n_h = n_h,
      total_n = sum(n_h),
      target_total_n = total_n,
      sample_size_gap = sum(n_h) - total_n, # how big the difference is between proposed and target sample size
      repaired_exactly = sum(n_h) == total_n,
      iterations = checked
    )
  }
  
  # Helper to score a solution
  score_result <- function(result) {
    c_change <- sum(abs(result$c_h - start_c_h))
    k_change <- sum(abs(result$k_h - start_k_h))
    
    c( # order matters - we want the n_total diff to be zero
      abs(result$sample_size_gap), # abs sample size gap
      c_change, # total change in ch
      k_change # total change in kh
    ) # we want to rather change kh than ch as the change in ch results in a more drastic change in n_total
  }
  
  # Compare two solutions
  # compare the scores between two solutions
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
  
  # Repair kh for fixed ch
  repair_k_for_c <- function(c_candidate, k_reference, checked) {
    # constrain candidate values
    c_candidate <- as.integer(pmax(lower_c, pmin(upper_c, c_candidate)))
    k_candidate <- as.integer(pmax(lower_k, pmin(upper_k, k_reference)))
    
    # best k rounded value thus far
    best_k <- k_candidate
    best_gap <- abs(sum(c_candidate * k_candidate) - total_n)
    
    step <- 0L
    
    while (step < max_iter) {
      step <- step + 1L
      
      # calculate the required adjustment
      current_total <- sum(c_candidate * k_candidate)
      diff <- total_n - current_total
      
      # save improvements
      if (abs(diff) < best_gap) {
        best_k <- k_candidate
        best_gap <- abs(diff)
      }
      
      # target is reached
      if (diff == 0L) {
        best_k <- k_candidate
        break
      }
      
      # current sample is too small
      if (diff > 0) {
        # find strata that can be increased without overshooting
        possible <- which(k_candidate < upper_k & c_candidate <= diff)
        
        if (length(possible) > 0) {
          j <- possible[which.max(c_candidate[possible])] # stratum with the largest ch that still fits in the gap
          jump <- min( # determines by how much we can increase kh
            floor(diff / c_candidate[j]), # largest number of kh increments that can be made without overshooting
            upper_k[j] - k_candidate[j],
            max(1L, as.integer(k_tolerance))
          )
          k_candidate[j] <- k_candidate[j] + as.integer(jump)
        } else { # fall back when every eligible increase will overshoot
          
          # increase kh for the smallest ch
          possible <- which(k_candidate < upper_k)
          if (length(possible) == 0) break
          
          j <- possible[which.min(c_candidate[possible])]
          trial_k <- k_candidate
          trial_k[j] <- trial_k[j] + 1L
          
          # only accepts the trial solution if it improves the gap
          if (abs(sum(c_candidate * trial_k) - total_n) >= best_gap) break
          k_candidate <- trial_k
        }
      } else { # current sample too large (identical to other case)
        diff_abs <- abs(diff)
        
        # find reductions that do not undershoot
        possible <- which(k_candidate > lower_k & c_candidate <= diff_abs)
        
        if (length(possible) > 0) {
          # prefer the largest reduction
          # reduce kh in the stratum with the largest ch that fits inside the excess
          j <- possible[which.max(c_candidate[possible])]
          
          jump <- min(
            floor(diff_abs / c_candidate[j]),
            k_candidate[j] - lower_k[j],
            max(1L, as.integer(k_tolerance))
          )
          
          k_candidate[j] <- k_candidate[j] - as.integer(jump)
        } else {
          # The fall back reduction
          # decrease kh in the stratum with the smallest ch
          possible <- which(k_candidate > lower_k)
          if (length(possible) == 0) break
          
          j <- possible[which.min(c_candidate[possible])]
          trial_k <- k_candidate
          trial_k[j] <- trial_k[j] - 1L
          
          if (abs(sum(c_candidate * trial_k) - total_n) >= best_gap) break
          k_candidate <- trial_k
        }
      }
    }
    
    # return the best kh for this ch
    make_result(c_candidate, best_k, checked)
  }
  
  c_candidates <- list(start_c_h)
  
  # generate nearby ch alternatives (in the search radius)
  # only one stratums ch is changed at a time
  for (j in seq_len(L)) {
    for (delta in seq.int(-c_search_radius, c_search_radius)) {
      if (delta == 0L) next
      
      candidate <- start_c_h
      candidate[j] <- candidate[j] + delta
      candidate <- as.integer(pmax(lower_c, pmin(upper_c, candidate)))
      c_candidates[[length(c_candidates) + 1L]] <- candidate
    }
  }
  
  # delete duplicate candidates
  c_keys <- vapply(c_candidates, paste, collapse = ":", character(1))
  c_candidates <- c_candidates[!duplicated(c_keys)]
  
  best <- NULL
  checked <- 0L
  
  # for each candidate
  for (candidate in c_candidates) {
    
    # holds the candidate ch fixed and repairs kh for that candidate
    checked <- checked + 1L
    repaired <- repair_k_for_c(candidate, start_k_h, checked)
    
    # compare the repaired allocation with the best so far 
    if (better_result(repaired, best)) {
      best <- repaired
    }
    
    # early stopping condition
    if (best$repaired_exactly && all(best$c_h == start_c_h)) {
      break
    }
    
    # max search limit
    if (checked >= max_iter) {
      break
    }
  }
  
  best$iterations <- checked
  best
}
