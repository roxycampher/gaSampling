repair_ckh <- function(c_h,
                       k_h,
                       total_n,
                       lower_c = NULL,
                       upper_c = NULL,
                       lower_k = NULL,
                       upper_k = NULL,
                       max_iter = 10000) {
  
  L <- length(c_h)
  
  # Defaults
  if (is.null(lower_c)) lower_c <- rep(1L, L) # every stratum must have at least one cluster
  if (is.null(upper_c)) upper_c <- rep(Inf, L) # need to set a true max
  if (is.null(lower_k)) lower_k <- rep(1L, L)
  if (is.null(upper_k)) upper_k <- rep(Inf, L)
  
  # make sure that the vectors all have the same length
  if (length(k_h) != L ||
      length(lower_c) != L ||
      length(upper_c) != L ||
      length(lower_k) != L ||
      length(upper_k) != L) {
    stop("c_h, k_h, and all bounds must have the same length.")
  }
  
  # round and enforce bounds
  c_h <- round(c_h)
  k_h <- round(k_h)
  
  c_h <- pmax(lower_c, pmin(upper_c, c_h)) # ensure the values are not too small
  k_h <- pmax(lower_k, pmin(upper_k, k_h))
  
  c_h <- as.integer(c_h)
  k_h <- as.integer(k_h)
  
  current_n <- sum(c_h * k_h)
  
  iter <- 0
  
  while (current_n != total_n && iter < max_iter) {
    iter <- iter + 1
    
    diff <- total_n - current_n
    
    if (diff > 0) {
      # Need to increase total sample size
      
      # Increasing c_h[j] by 1 increases n by k_h[j]
      possible_c <- which(c_h < upper_c & k_h <= diff)
      
      # Increasing k_h[j] by 1 increases n by c_h[j]
      possible_k <- which(k_h < upper_k & c_h <= diff)
      
      if (length(possible_c) == 0 && length(possible_k) == 0) break
      
      # Table of possible moves
      moves <- data.frame(
        type = character(0),
        j = integer(0),
        increase = numeric(0)
      )
      
      if (length(possible_c) > 0) { # add the possible moves for c
        moves <- rbind(
          moves,
          data.frame(
            type = "c",
            j = possible_c,
            increase = k_h[possible_c]
          )
        )
      }
      
      if (length(possible_k) > 0) { # add the possible moves for k 
        moves <- rbind(
          moves,
          data.frame(
            type = "k",
            j = possible_k,
            increase = c_h[possible_k]
          )
        )
      }
      
      # Choose the largest increase that does not overshoot n_total (Greedy)
      best_move <- moves[which.max(moves$increase), ]
      
      if (best_move$type == "c") {
        c_h[best_move$j] <- c_h[best_move$j] + 1L
      } else {
        k_h[best_move$j] <- k_h[best_move$j] + 1L
      }
      
    } else {
      # Need to decrease total sample size
      
      diff_abs <- abs(diff)
      
      # Decreasing c_h[j] by 1 decreases n by k_h[j]
      possible_c <- which(c_h > lower_c & k_h <= diff_abs)
      
      # Decreasing k_h[j] by 1 decreases n by c_h[j]
      possible_k <- which(k_h > lower_k & c_h <= diff_abs)
      
      if (length(possible_c) == 0 && length(possible_k) == 0) break
      
      moves <- data.frame(
        type = character(0),
        j = integer(0),
        decrease = numeric(0)
      )
      
      if (length(possible_c) > 0) {
        moves <- rbind(
          moves,
          data.frame(
            type = "c",
            j = possible_c,
            decrease = k_h[possible_c]
          )
        )
      }
      
      if (length(possible_k) > 0) {
        moves <- rbind(
          moves,
          data.frame(
            type = "k",
            j = possible_k,
            decrease = c_h[possible_k]
          )
        )
      }
      
      # Choose the largest decrease that does not undershoot total_n
      best_move <- moves[which.max(moves$decrease), ]
      
      if (best_move$type == "c") {
        c_h[best_move$j] <- c_h[best_move$j] - 1L
      } else {
        k_h[best_move$j] <- k_h[best_move$j] - 1L
      }
    }
    
    current_n <- sum(c_h * k_h)
  }
  
  list(
    c_h = as.integer(c_h),
    k_h = as.integer(k_h),
    n_h = as.integer(c_h * k_h),
    total_n = sum(c_h * k_h),
    repaired_exactly = sum(c_h * k_h) == total_n,
    iterations = iter
  )
}