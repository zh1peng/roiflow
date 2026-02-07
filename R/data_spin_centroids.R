#' Left-Hemisphere Centroid Coordinates for Spin Tests
#'
#' @description
#' Centroid coordinates for 34 left-hemisphere Desikan cortical regions,
#' used to generate spherical spin permutations in \code{\link{spin_correlation}}.
#'
#' Columns:
#' \itemize{
#'   \item \code{Row}: ROI label (e.g., \code{"L_superiorfrontal"})
#'   \item \code{centroid1}, \code{centroid2}, \code{centroid3}: 3D coordinates
#' }
#'
#' @format A data frame with 34 rows and 4 variables.
#' @usage data(lh_centroid_data)
"lh_centroid_data"

#' Right-Hemisphere Centroid Coordinates for Spin Tests
#'
#' @description
#' Centroid coordinates for 34 right-hemisphere Desikan cortical regions,
#' used to generate spherical spin permutations in \code{\link{spin_correlation}}.
#'
#' Columns:
#' \itemize{
#'   \item \code{Row}: ROI label (e.g., \code{"R_superiorfrontal"})
#'   \item \code{centroid1}, \code{centroid2}, \code{centroid3}: 3D coordinates
#' }
#'
#' @format A data frame with 34 rows and 4 variables.
#' @usage data(rh_centroid_data)
"rh_centroid_data"
