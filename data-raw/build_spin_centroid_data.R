# Build packaged centroid datasets for spin tests.
# Run from project root: Rscript data-raw/build_spin_centroid_data.R

lh_centroid_data <- utils::read.csv('data-raw/lh_centroid.csv', stringsAsFactors = FALSE)
rh_centroid_data <- utils::read.csv('data-raw/rh_centroid.csv', stringsAsFactors = FALSE)

# Keep canonical column names expected by spin utilities.
if ('row' %in% names(lh_centroid_data) && !('Row' %in% names(lh_centroid_data))) {
  names(lh_centroid_data)[names(lh_centroid_data) == 'row'] <- 'Row'
}
if ('row' %in% names(rh_centroid_data) && !('Row' %in% names(rh_centroid_data))) {
  names(rh_centroid_data)[names(rh_centroid_data) == 'row'] <- 'Row'
}

required_cols <- c('Row', 'centroid1', 'centroid2', 'centroid3')
missing_lh <- setdiff(required_cols, names(lh_centroid_data))
missing_rh <- setdiff(required_cols, names(rh_centroid_data))
if (length(missing_lh) > 0) stop(sprintf('lh_centroid_data missing columns: %s', paste(missing_lh, collapse = ', ')))
if (length(missing_rh) > 0) stop(sprintf('rh_centroid_data missing columns: %s', paste(missing_rh, collapse = ', ')))

lh_centroid_data <- lh_centroid_data[, required_cols]
rh_centroid_data <- rh_centroid_data[, required_cols]

save(lh_centroid_data, file = 'data/lh_centroid_data.rda', compress = 'bzip2')
save(rh_centroid_data, file = 'data/rh_centroid_data.rda', compress = 'bzip2')
if (file.exists('data/spin_centroids.rda')) file.remove('data/spin_centroids.rda')
message('Wrote data/lh_centroid_data.rda and data/rh_centroid_data.rda')
