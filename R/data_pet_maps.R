#' PET Maps in Desikan-Killiany 68 Space
#'
#' @description
#' Group-level PET receptor/transporter maps in DK-68 ROI space.
#' Rows correspond to cortical ROIs (left then right hemisphere), and columns
#' correspond to PET targets (e.g., \code{D1}, \code{5-HTT}, \code{DAT}).
#'
#' These maps are derived from \code{neuromaps} annotations using the
#' \code{data-raw/prepare_PET_maps.py} pipeline, then summarized as weighted
#' means across studies per target.
#'
#' @format A numeric matrix with 68 rows (ROIs) and 19 columns (PET maps).
#' Row names are ROI labels (e.g., \code{"L_bankssts"}, \code{"R_insula"}).
#' @usage data(pet_maps_dk68)
"pet_maps_dk68"

#' Metadata for Packaged PET Maps
#'
#' @description
#' Map-level metadata for \code{pet_maps_dk68}, including target label,
#' modality, and processing/parcellation assumptions.
#'
#' \itemize{
#'   \item \code{pet_map}: map short name used in \code{pet_select}
#'   \item \code{target}: human-readable PET target label
#'   \item \code{modality}: \code{"receptor"} or \code{"transporter"}
#'   \item \code{tracer}: tracer string if available in source metadata (currently NA)
#'   \item \code{source}: data provenance
#'   \item \code{space}: surface/space conversion summary
#'   \item \code{parcellation}: ROI atlas assumption
#'   \item \code{preprocessing}: processing summary from data-raw pipeline
#' }
#'
#' @format A data frame with 19 rows and 8 variables.
#' @usage data(pet_maps_meta)
"pet_maps_meta"
