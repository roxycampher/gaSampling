gcd_two <- function(a, b) {
  a <- abs(a)
  b <- abs(b)
  
  if (b == 0) return(a)
  
  Recall(b, a %% b) # recursive algorithm to calcuate the GCD
}

gcd_multiple <- function(x) {
  x <- abs(x)
  Reduce(gcd_two, x)
}

# Checks whether the requested n_total is feasible for the given c_h
check_total_n_feasible <- function(total_n, c_h, lower = NULL, upper = NULL) {
  if (is.null(lower)) lower <- rep(1L, length(c_h))
  if (is.null(upper)) upper <- rep(Inf, length(c_h))
  
  g <- gcd_multiple(c_h)
  
  gcd_ok <- total_n %% g == 0 # check if n_total is divisible by gcd of c_h
  range_ok <- total_n >= sum(c_h * lower) && total_n <= sum(c_h * upper) # check if n_total is within the possible range
  
  feasible <- gcd_ok && range_ok # both checks must pass
  
  list(
    feasible = feasible,
    gcd_ok = gcd_ok,
    range_ok = range_ok,
    gcd_c_h = g,
    min_n = sum(c_h * lower), # smallest possible sample size under the given c_h
    max_n = sum(c_h * upper), # largest possible sample size under the given c_h
    total_n = total_n,
    message = if (feasible) {
      "Exact total_n is feasible for this c_h allocation."
    } else {
      paste0(
        "Exact total_n may be impossible. total_n = ", total_n,
        ", gcd(c_h) = ", g,
        ", min possible n = ", sum(c_h * lower),
        ", max possible n = ", sum(c_h * upper), "."
      )
    }
  )
}