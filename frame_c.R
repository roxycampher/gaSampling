#' Affordable cluster count function
#'
#' Calculate the total number of available clusters affordable within a budget
#'
#' Count distinct stratum-cluster pairs in the survey data and calculate
#' how many clusters can be afforded after reserving the cost of sampling
#' the specified total number of secondary sampling units.
#' Return the smaller of the available and affordable cluster counts.
#'
#' @param survey_data A data frame containing columns named strata
#'   and cluster, which identify the stratum and primary sampling unit
#'   for each observation.
#' @param max_cost A single finite non-negative numeric value giving
#'   the total available budget.
#' @param unit_cost A single finite non-negative numeric value giving
#'   the cost per sampled secondary sampling unit.
#' @param cluster_cost A single finite positive numeric value giving
#'   the cost per sampled primary sampling unit (cluster).
#' @param total_n A single non-negative integer giving the total number
#'   of secondary sampling units to be sampled.
#'
#' @return A single numeric value giving the smaller of the number of
#'   available clusters and the number affordable within the remaining
#'   budget. May be negative if unit sampling costs exceed the budget.
#'
#' @details
#'   Available clusters are counted as distinct combinations of strata
#'   and cluster, allowing cluster identifiers to repeat across strata.
#'   These identifiers should not contain missing values.
#'
#'   The affordable cluster count is calculated as
#'   floor((max_cost - unit_cost * total_n) / cluster_cost).
#'   Unit and cluster costs are assumed to be constant across strata.
#'
#'   The result is a total across all strata, rather than a separate
#'   count for each stratum. It does not check whether the available
#'   clusters contain enough units to achieve total_n.
#'   The function does not validate inputs or restrict negative
#'   results to zero.
#'
#' @examples
#' survey_data <- data.frame(
#'   strata = c(1, 1, 1, 2, 2, 2),
#'   cluster = c(1, 1, 2, 1, 2, 3)
#' )
#'
#' frame_c(
#'   survey_data = survey_data,
#'   max_cost = 1000,
#'   unit_cost = 5,
#'   cluster_cost = 100,
#'   total_n = 120
#' )
frame_c <- function(survey_data, max_cost, unit_cost, cluster_cost, total_n) {
  
  available_c <- nrow(unique(survey_data[c("strata", "cluster")]))
  affordable_c <- floor((max_cost - unit_cost * total_n) / cluster_cost)
  
  return(min(available_c, affordable_c))
}
