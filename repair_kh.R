repair_kh <- function(x, c_h, total_n, lower = NULL, upper = NULL) {
  # round to integers
  k_h <- round(x)
  
  # apply lower and upper bounds (defaults)
  if (is.null(lower)) lower <- rep(1L, length(k_h))
  if (is.null(upper)) upper <- rep(Inf, length(k_h))
  
  k_h <- pmax(lower, pmin(upper, k_h)) # clamps k_h values between its lower and upper bounds
  
  # current achieved total sample size
  current_n <- sum(c_h * k_h)
  
  # repair until total_n is matched or no move possible
  max_iter <- 10000
  iter <- 0
  
  while (current_n != total_n && iter < max_iter) {
    iter <- iter + 1
    
    diff <- total_n - current_n
    
    if (diff > 0) {
      # Need to increase sample size
      possible <- which(k_h < upper & c_h <= diff) # can we find a stratum where we can increase k_h by 1
      
      if (length(possible) == 0) break # if no valid stratum can be used the loop ends
      
      # Choose stratum with the largest ch
      j <- possible[which.max(c_h[possible])]
      
      k_h[j] <- k_h[j] + 1
      
    } else {
      # Need to decrease sample size
      possible <- which(k_h > lower & c_h <= abs(diff)) # find strata where we can decrease k_h by 1
      
      if (length(possible) == 0) break # no valid decrease is possible the loop stops
      
      # Choose stratum that gives largest feasible decrease
      j <- possible[which.max(c_h[possible])]
      
      k_h[j] <- k_h[j] - 1
    }
    
    current_n <- sum(c_h * k_h) # update the current total
  }
  
  return(as.integer(k_h)) # returns the repaired k_h vector as integers
}