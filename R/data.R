#' ECD census extract
#'
#' A small extract of the 2021 South African early childhood development (ECD)
#' census, used to illustrate the allocation workflow. It contains all ECD
#' centres in a random sample of up to 40 wards per province, with at least two
#' centres in each ward. It is intended for illustration only and should not be
#' used for inference about the full census.
#'
#' @format A data frame with one row per ECD centre and the following columns:
#' \describe{
#'   \item{id}{ECD centre identifier (secondary sampling unit).}
#'   \item{province}{Province name (stratum).}
#'   \item{id_prov}{Province code.}
#'   \item{id_dc}{District council code.}
#'   \item{dc}{District council name.}
#'   \item{id_mn}{Municipality code.}
#'   \item{mn}{Municipality name.}
#'   \item{id_ward}{Ward identifier (cluster, the primary sampling unit).}
#'   \item{count_children_present}{Number of children present at the centre
#'     (response variable).}
#'   \item{capacity}{Capacity of the centre.}
#'   \item{count_staff_all}{Total number of staff at the centre.}
#'   \item{fees_amount}{Fees charged by the centre.}
#' }
#' @source Department of Basic Education, South Africa. Early Childhood
#'   Development Census 2021 (reference ID zaf-dbe-ecdc-2021-v2), public-use
#'   version, accessed through DataFirst, University of Cape Town.
"ecd_subset"