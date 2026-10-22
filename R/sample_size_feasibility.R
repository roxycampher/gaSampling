#' Greatest common divisor function
#'
#' Calculate the greatest common divisor of two integers
#'
#' Use the Euclidean algorithm to find the largest positive integer
#' that divides both inputs without a remainder.
#'
#' @param a A single finite integer-valued number.
#' @param b A single finite integer-valued number.
#'
#' @return A single numeric value giving the non-negative greatest
#'   common divisor of a and b. Returns zero when both inputs are zero.
#'
#' @details
#'   Negative inputs are converted to their absolute values.
#'   The calculation uses recursion and does not validate inputs.
#'
#' @examples
#' gcd_two(a = 12, b = 18)
#' gcd_two(a = -12, b = 0)
#' @keywords internal
#' @noRd
gcd_two <- function(a, b) {
  a <- abs(a)
  b <- abs(b)
  
  if (b == 0) return(a)
  
  Recall(b, a %% b) # recursive algorithm to calcuate the GCD
}

#' Multiple-value greatest common divisor function
#'
#' Calculate the greatest common divisor of a vector of integers
#'
#' Apply gcd_two repeatedly to find the greatest common divisor
#' shared by all values in the supplied vector.
#'
#' @param x A non-empty numeric vector of finite integer-valued numbers.
#'
#' @return A single numeric value giving the non-negative greatest
#'   common divisor of the supplied values. Returns zero when all
#'   values are zero.
#'
#' @details
#'   Negative inputs are converted to their absolute values.
#'   A vector containing one value returns the absolute value of that
#'   element. The function requires gcd_two to be available and does
#'   not validate inputs.
#'
#' @examples
#' gcd_multiple(x = c(12, 18, 24))
#' gcd_multiple(x = c(-10, 20, 30))
#' @keywords internal
#' @noRd
gcd_multiple <- function(x) {
  x <- abs(x)
  Reduce(gcd_two, x)
}

#' Minimum sample size feasibility function
#'
#' Check whether a minimum sample size can be reached
#'
#' Determine whether fixed numbers of sampled clusters and integer
#' bounds on the number of sampled units per cluster allow a total
#' sample size that meets or exceeds the required minimum.
#'
#' @param total_n A single finite non-negative integer-valued number
#'   giving the required minimum total sample size.
#' @param c_h Numeric vector giving the fixed number of primary sampling
#'   units (clusters) sampled in each stratum. Values should be positive
#'   integers.
#' @param lower Numeric vector giving the lower bound on the number of
#'   secondary sampling units sampled per cluster in each stratum.
#'   Values should be finite and positive. Defaults to one per stratum
#'   when NULL.
#' @param upper Numeric vector giving the upper bound on the number of
#'   secondary sampling units sampled per cluster in each stratum.
#'   Values should be positive; Inf indicates no upper limit.
#'   Defaults to Inf per stratum when NULL.
#'
#' @return A named list containing:
#'   \item{feasible}{Logical value indicating whether the bounds are
#'     consistent and the required minimum does not exceed the maximum
#'     attainable sample size.}
#'   \item{bounds_ok}{Logical value indicating whether each rounded
#'     lower bound is no greater than its rounded upper bound.}
#'   \item{range_ok}{Logical value indicating whether total_n is no
#'     greater than max_n.}
#'   \item{min_n}{Minimum total sample size implied by the rounded
#'     lower bounds.}
#'   \item{max_n}{Maximum total sample size implied by the rounded
#'     upper bounds. May be Inf.}
#'   \item{total_n}{The supplied minimum sample size requirement.}
#'   \item{message}{A character string explaining the feasibility result.}
#'
#' @details
#'   The vectors c_h, lower and upper should have the same length,
#'   with entries ordered consistently by stratum.
#'   Lower bounds are rounded up and upper bounds are rounded down
#'   to enforce integer numbers of sampled units per cluster.
#'   Total sample size is calculated as sum(c_h * k_h).
#'
#'   This function checks whether the minimum requirement can be met
#'   or exceeded, rather than whether an exact total can be attained.
#'   A requirement below min_n is therefore feasible when the bounds
#'   are consistent. The function does not construct an allocation
#'   or validate input types, lengths or positivity.
#'
#' @examples
#' check_total_n_feasible(
#'   total_n = 150,
#'   c_h = c(10, 15),
#'   lower = c(2, 3),
#'   upper = c(5, 8)
#' )
#' @export
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