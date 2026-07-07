prepare_strata <- function(data, # might have to change for the joint
                           strata_var,
                           cluster_var,
                           unit_var,
                           y_var,
                           cost_cluster = 1,
                           cost_unit = 1) {
  
  # Keep only the needed variables
  dat <- data[, c(strata_var, cluster_var, unit_var, y_var)]
  
  names(dat) <- c("stratum", "cluster", "unit", "y")
  
  # Remove missing values
  dat <- dat[!is.na(dat$stratum) &
               !is.na(dat$cluster) &
               !is.na(dat$unit) &
               !is.na(dat$y), ]
  
  # Total number of units
  N_total <- length(unique(dat$unit))
  
  # Cluster-level summaries
  cluster_summary <- dat |>
    dplyr::group_by(stratum, cluster) |>
    dplyr::summarise(
      cluster_mean = mean(y, na.rm = TRUE),
      cluster_var = stats::var(y, na.rm = TRUE),
      cluster_n = dplyr::n(),
      .groups = "drop"
    )
  
  # Stratum-level summaries
  strata <- dat |>
    dplyr::group_by(stratum) |>
    dplyr::summarise(
      N_h = dplyr::n_distinct(unit),       # number of units in stratum h
      C_h = dplyr::n_distinct(cluster),    # number of clusters available in stratum h
      .groups = "drop"
    )
  
  # Between-cluster variance
  between_var <- cluster_summary |>
    dplyr::group_by(stratum) |>
    dplyr::summarise(
      sigma2_bh = stats::var(cluster_mean, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Within-cluster variance
  within_var <- cluster_summary |>
    dplyr::group_by(stratum) |>
    dplyr::summarise(
      sigma2_wh = weighted.mean(
        x = cluster_var,
        w = cluster_n,
        na.rm = TRUE
      ),
      .groups = "drop"
    )
  
  # Join everything together
  strata <- strata |>
    dplyr::left_join(between_var, by = "stratum") |>
    dplyr::left_join(within_var, by = "stratum") |>
    dplyr::mutate(
      W_h = N_h / N_total,
      cost_cluster = cost_cluster,
      cost_unit = cost_unit
    )
  
  # Replace NA variance values if a stratum has too few clusters
  strata$sigma2_bh[is.na(strata$sigma2_bh)] <- 0
  strata$sigma2_wh[is.na(strata$sigma2_wh)] <- 0
  
  strata
}