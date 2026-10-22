#' Two-stage variance function
#' 
#' Calculate variance for a stratified two-stage allocation
#'
#' Calculate the variance of the overall mean by summing weighted
#' between-cluster and within-cluster variance contributions across strata.
#'
#' @param c_h Numeric vector giving the number of primary sampling units
#'   (clusters) sampled in each stratum. Values must be positive.
#' @param k_h Numeric vector giving the number of secondary sampling units
#'   sampled per cluster in each stratum. Values must be positive.
#' @param W_h Numeric vector of stratum population weights, calculated as
#'   N_h / sum(N_h), where N_h is the population size of stratum h.
#' @param sigma2_bh Numeric vector of between-cluster variance components
#'   for each stratum. Values must be non-negative.
#' @param sigma2_wh Numeric vector of within-cluster variance components
#'   for each stratum. Values must be non-negative.
#'
#' @return A single numeric value giving the variance for the supplied
#'   allocation.
#'
#' @details
#'   All inputs should be finite vectors of the same length, with entries
#'   ordered consistently by stratum. The function does not validate inputs.
#'   The calculation does not include finite population corrections.
#'
#' @examples
#' variance_twostage(
#'   c_h = c(10, 15),
#'   k_h = c(5, 8),
#'   W_h = c(0.4, 0.6),
#'   sigma2_bh = c(2, 3),
#'   sigma2_wh = c(6, 9)
#' )
#' @export

variance_twostage <- function(c_h, k_h, W_h, sigma2_bh, sigma2_wh) {
  sum(W_h^2 * ((sigma2_bh / c_h) + (sigma2_wh / (c_h * k_h))))
}