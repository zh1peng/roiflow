# Build packaged PET datasets for DK-68 analyses.
# Run from project root: Rscript data-raw/build_pet_data.R

pet_csv <- "data-raw/PET_DK68_weighted_mean.csv"
if (!file.exists(pet_csv)) {
  stop(sprintf("PET source file not found: %s", pet_csv), call. = FALSE)
}

pet_tbl <- utils::read.csv(pet_csv, stringsAsFactors = FALSE, check.names = FALSE)
required_cols <- c("RegionName")
missing_cols <- setdiff(required_cols, names(pet_tbl))
if (length(missing_cols) > 0) {
  stop(
    sprintf("PET source is missing required columns: %s", paste(missing_cols, collapse = ", ")),
    call. = FALSE
  )
}

map_names <- setdiff(names(pet_tbl), "RegionName")
if (length(map_names) == 0) {
  stop("No PET map columns found in PET_DK68_weighted_mean.csv.", call. = FALSE)
}

pet_maps_dk68 <- as.matrix(pet_tbl[, map_names, drop = FALSE])
storage.mode(pet_maps_dk68) <- "double"
rownames(pet_maps_dk68) <- pet_tbl$RegionName
colnames(pet_maps_dk68) <- map_names

target_lookup <- c(
  VAChT = "Vesicular acetylcholine transporter",
  DAT = "Dopamine transporter",
  NET = "Norepinephrine transporter",
  `5-HTT` = "Serotonin transporter",
  D1 = "Dopamine D1 receptor",
  D2 = "Dopamine D2 receptor",
  `5-HT1a` = "Serotonin 1A receptor",
  `5-HT1b` = "Serotonin 1B receptor",
  `5-HT2a` = "Serotonin 2A receptor",
  `5-HT4` = "Serotonin 4 receptor",
  `5-HT6` = "Serotonin 6 receptor",
  mGluR5 = "Metabotropic glutamate receptor 5",
  GABAa = "GABA-A receptor",
  H3 = "Histamine H3 receptor",
  a4b2 = "Nicotinic acetylcholine receptor alpha4beta2",
  MOR = "Mu-opioid receptor",
  KOR = "Kappa-opioid receptor",
  CB1 = "Cannabinoid receptor type 1",
  M1 = "Muscarinic acetylcholine receptor M1"
)

transporter_maps <- c("VAChT", "DAT", "NET", "5-HTT")
pet_maps_meta <- data.frame(
  pet_map = map_names,
  target = unname(target_lookup[map_names]),
  modality = ifelse(map_names %in% transporter_maps, "transporter", "receptor"),
  tracer = NA_character_,
  source = "Derived from neuromaps annotations via data-raw/prepare_PET_maps.py",
  space = "fsaverage10k to DK68",
  parcellation = "Desikan-Killiany 68",
  preprocessing = "Negative values clipped at vertex-level; DK values min-max scaled; weighted mean by study total_n",
  stringsAsFactors = FALSE
)

if (anyNA(pet_maps_meta$target)) {
  missing_targets <- pet_maps_meta$pet_map[is.na(pet_maps_meta$target)]
  stop(
    sprintf("Missing metadata target label for map(s): %s", paste(missing_targets, collapse = ", ")),
    call. = FALSE
  )
}

save(pet_maps_dk68, file = "data/pet_maps_dk68.rda", compress = "bzip2")
save(pet_maps_meta, file = "data/pet_maps_meta.rda", compress = "bzip2")
message("Wrote data/pet_maps_dk68.rda and data/pet_maps_meta.rda")
