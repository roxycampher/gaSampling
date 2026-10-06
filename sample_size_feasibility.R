#' Calculate the greatest common divisor (GCD) of two integers
#' using the recursive Euclidean algorithm.
#'
#' Inputs:
#'   a, b: Finite integer-valued numeric scalars. Negative values are converted to their absolute values.
#'
#' Output:
#'   A single non-negative numeric value giving the GCD.
#'   Returns 0 if both inputs are 0.
gcd_two <- function(a, b) {
  a <- abs(a)
  b <- abs(b)
  
  if (b == 0) return(a)
  
  Recall(b, a %% b) # recursive algorithm to calcuate the GCD
}

#' Calculate the greatest common divisor (GCD) of a vector of integers
#' by repeatedly applying gcd_two().
#'
#' Inputs:
#'   x: A non-empty numeric vector of finite integer-valued numbers.
#'      Negative values are converted to their absolute values.
#'
#' Output:
#'   A single non-negative numeric value giving the GCD of all entries.
#'   Returns 0 if all entries are 0.
#'
#' Dependencies:
#'   Requires gcd_two() to be defined.
gcd_multiple <- function(x) {
  x <- abs(x)
  Reduce(gcd_two, x)
}

# Check whether the minimum sample size is attainable
# for fixed c_h, subject to integer k_h bounds.
check_total_n_feasible <- function(total_n, c_h,
                                   lower = NULL, upper = NULL) {
  if (is.null(lower)) lower <- rep(1L, length(c_h))
  if (is.null(upper)) upper <- rep(Inf, length(c_h))
  
  lower <- ceiling(lower)
  upper <- floor(upper)
  
  bounds_ok <- all(lower <= upper)
  min_n <- sum(c_h * lower)
  max_n <- sum(c_h * upper)
  
  range_ok <- total_n <= max_n
  feasible <- bounds_ok && range_ok
  
  list(
    feasible = feasible,
    bounds_ok = bounds_ok,
    range_ok = range_ok,
    min_n = min_n,
    max_n = max_n,
    total_n = total_n,
    message = if (!bounds_ok) {
      "The lower and upper bounds on k_h are inconsistent."
    } else if (feasible) {
      "The minimum sample size is attainable under these allocation bounds."
    } else {
      "The minimum sample size exceeds the maximum attainable sample size."
    }
  )
}