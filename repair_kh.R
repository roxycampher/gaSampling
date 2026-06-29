repair_kh <- function(x, c_h, total_n, lower = NULL, upper = NULL) {
  # round to integers
  k_h <- round(x)
  
  # apply lower and upper bounds
  if (is.null(lower)) lower <- rep(1L, length(k_h))
  if (is.null(upper)) upper <- rep(Inf, length(k_h))
  
  k_h <- pmax(lower, pmin(upper, k_h))
  
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
      possible <- which(k_h < upper & c_h <= diff)
      
      if (length(possible) == 0) break
      
      # Choose stratum that gives largest feasible increase
      j <- possible[which.max(c_h[possible])]
      
      k_h[j] <- k_h[j] + 1
      
    } else {
      # Need to decrease sample size
      possible <- which(k_h > lower & c_h <= abs(diff))
      
      if (length(possible) == 0) break
      
      # Choose stratum that gives largest feasible decrease
      j <- possible[which.max(c_h[possible])]
      
      k_h[j] <- k_h[j] - 1
    }
    
    current_n <- sum(c_h * k_h)
  }
  
  return(as.integer(k_h))
}