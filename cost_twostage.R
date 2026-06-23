cost_twostage <- function(c_h, k_h, cost_cluster, cost_unit) {
  sum(cost_cluster * c_h + cost_unit * c_h * k_h)
}