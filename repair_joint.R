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
  
}