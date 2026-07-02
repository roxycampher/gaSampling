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
      
    } else {
      # Need to decrease total sample size
      
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