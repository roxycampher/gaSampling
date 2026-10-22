#' Two-stage cost function
#'
#' Calculate the total cost of a stratified two-stage allocation
#'
#' Sum the costs of sampling primary sampling units (clusters) and
#' secondary sampling units within those clusters across all strata.
#'
#' @param c_h Numeric vector giving the number of primary sampling units
#'   (clusters) sampled in each stratum. Values should be non-negative
#'   integers.
#' @param k_h Numeric vector giving the number of secondary sampling
#'   units sampled per cluster in each stratum. Values should be
#'   non-negative integers.
#' @param cost_cluster Numeric vector giving the cost per sampled
#'   cluster in each stratum, or a single value applied to all strata.
#'   Values must be finite and non-negative.
#' @param cost_unit Numeric vector giving the cost per sampled secondary
#'   unit in each stratum, or a single value applied to all strata.
#'   Values must be finite and non-negative.
#'
#' @return A single numeric value giving the total sampling cost for
#'   the supplied allocation, in the same currency units as the inputs.
#'
#' @details
#'   Allocation vectors should be finite and have the same length,
#'   with entries ordered consistently by stratum. Cost vectors should
#'   have one entry per stratum unless a single common cost is supplied.
#'
#'   Each stratum contributes cost_cluster * c_h in cluster costs and
#'   cost_unit * c_h * k_h in unit costs. Cluster costs are assumed to
#'   be additional to unit costs.
#'   The function does not validate inputs or check budget constraints.
#'
#' @examples
#' cost_twostage(
#'   c_h = c(10, 15),
#'   k_h = c(5, 8),
#'   cost_cluster = c(100, 120),
#'   cost_unit = c(5, 6)
#' )
#' @export
cost_twostage <- function(c_h, k_h, cost_cluster, cost_unit) {
  sum(cost_cluster * c_h + cost_unit * c_h * k_h)
}