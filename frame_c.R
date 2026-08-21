#-------------------------------------------------------------------------------
# Function to calculate the number of clusters that is realistic per stratum
#-------------------------------------------------------------------------------

frame_c <- function(survey_data, max_cost, unit_cost, cluster_cost, total_n) {
  
  available_c <- nrow(unique(survey_data[c("strata", "cluster")]))
  affordable_c <- floor((max_cost - unit_cost * total_n) / cluster_cost)
  
  return(min(available_c, affordable_c))
}
