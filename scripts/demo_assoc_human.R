#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("devtools is required to run this demo from source. Install with install.packages('devtools').", call. = FALSE)
  }
  devtools::load_all(".")
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for this demo. Install with install.packages('ggplot2').", call. = FALSE)
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("jsonlite is required for this demo. Install with install.packages('jsonlite').", call. = FALSE)
  }
})

# Human demo (association):
# End-to-end (prep -> ComBat -> association scan -> tables -> plots -> manifest).
#
# Goal: association scan across ROIs:
#   ROI ~ Age + Sex

# ------------------------------ Parameters -----------------------------------
input_file <- file.path("data", "sample.csv")
out_root <- file.path("outputs", "assoc_human")

# Columns (hard-coded for the sample dataset)
site_col <- "Site"
age_col <- "Age"
sex_col <- "Sex"
icv_col <- "ICV"

# ROI detection
roi_regex <- "^(L_|R_)"
roi_drop_suffix <- "_thickavg"

# Example ROI for scatter plot
scatter_roi <- "L_superiorfrontal_thickavg"
scatter_value <- "marginal" # "raw" or "marginal"

if (!file.exists(input_file)) stop(sprintf("Input file not found: %s", input_file), call. = FALSE)

dir.create(file.path(out_root, "data"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_root, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_root, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_root, "logs"), showWarnings = FALSE, recursive = TRUE)

timings <- list()
t0_all <- Sys.time()

## Step 0: Load Data -----------------------------------------------------------
cat("\n## Step 0: Load Data\n")
t0 <- Sys.time()
df0 <- read.csv(input_file, stringsAsFactors = FALSE)
t1 <- Sys.time()
timings$load <- list(start = as.character(t0), end = as.character(t1),
                     seconds = as.numeric(difftime(t1, t0, units = "secs")))

required <- c(site_col, age_col, sex_col, icv_col)
missing <- setdiff(required, names(df0))
if (length(missing) > 0) stop(sprintf("Missing required columns: %s", paste(missing, collapse = ", ")), call. = FALSE)

roi_cols <- grep(roi_regex, names(df0), value = TRUE)
if (length(roi_cols) == 0) stop(sprintf("No ROI columns found with roi_regex='%s'.", roi_regex), call. = FALSE)

## Step 1: Prep ----------------------------------------------------------------
cat("\n## Step 1: Prep\n")
t0 <- Sys.time()

df0[[site_col]] <- factor(as.character(df0[[site_col]]))
df0[[sex_col]] <- factor(as.character(df0[[sex_col]]))

sex_ref <- if ("0" %in% levels(df0[[sex_col]])) "0" else levels(df0[[sex_col]])[1]
factor_levels <- list()
factor_levels[[sex_col]] <- sex_ref

prep_spec_obj <- prep_spec(
  roi_regex = roi_regex,
  site_col = site_col,
  factor_levels = factor_levels,
  na_action = "drop_required",
  outlier_action = "none"
)

prep_res <- prep(df0, spec = prep_spec_obj, return = "result", verbose = FALSE)
df_prep <- prep_res$data

t1 <- Sys.time()
timings$prep <- list(start = as.character(t0), end = as.character(t1),
                     seconds = as.numeric(difftime(t1, t0, units = "secs")))

## Step 2: ComBat --------------------------------------------------------------
cat("\n## Step 2: ComBat\n")
t0 <- Sys.time()

combat_spec_obj <- combat_spec(
  roi_cols = roi_cols,
  roi_regex = roi_regex,
  batch_col = site_col,
  mod_formula = stats::as.formula(paste0("~ ", age_col, " + ", sex_col, " + ", icv_col)),
  missing_fit = "impute",
  single_batch = "error",
  keep_inputs = TRUE
)

combat_fit_obj <- combat_fit_df(df_prep, spec = combat_spec_obj, verbose = TRUE)
combat_apply_res <- combat_apply_df(combat_fit_obj, df_prep, verbose = TRUE)
df_combat <- combat_apply_res$data

write.csv(df_combat, file.path(out_root, "data", "sample_combat.csv"), row.names = FALSE)

t1 <- Sys.time()
timings$combat <- list(start = as.character(t0), end = as.character(t1),
                       seconds = as.numeric(difftime(t1, t0, units = "secs")))

## Step 3: Association Scan (Age + Sex) ---------------------------------------
cat("\n## Step 3: Association Scan\n")
t0 <- Sys.time()

analysis_spec_obj <- analysis_spec(
  covariates = c(sex_col),
  model_engine = "lm",
  diagnostics = "light",
  effect_size = "standardized_beta",
  effect_size_from = "model",
  p_adjust = list(method = "fdr", scope = "within_call",
                  family_id = "ct_age_assoc", family_desc = "CT ~ Age + Sex"),
  ref_levels = setNames(list(sex_ref), sex_col)
)

res_assoc <- assoc_scan(df_combat, y_vars = roi_cols, x = age_col, spec = analysis_spec_obj, return = "result")
res_tbl <- res_assoc$data

write.csv(res_tbl, file.path(out_root, "tables", "ct_age_assoc_raw.csv"), row.names = FALSE)
if (!is.null(res_assoc$diagnostics)) {
  write.csv(res_assoc$diagnostics, file.path(out_root, "tables", "ct_age_assoc_diagnostics.csv"), row.names = FALSE)
}

t1 <- Sys.time()
timings$assoc_scan <- list(start = as.character(t0), end = as.character(t1),
                           seconds = as.numeric(difftime(t1, t0, units = "secs")))

## Step 4: Plots ---------------------------------------------------------------
cat("\n## Step 4: Plots\n")
t0 <- Sys.time()

qc_plot <- plot_qc_pvalue_hist(
  res_tbl,
  p_col = "p_adj",
  title = "QC: Age Association (FDR p-value distribution)",
  plot_spec = list(theme = "minimal")
)
ggplot2::ggsave(file.path(out_root, "figures", "qc_pvalue_hist.png"), plot = qc_plot, width = 6, height = 4, bg = "white")

top_rois_plot <- plot_top_rois(
  res_tbl,
  metric = "estimate",
  n = 20,
  label_col = "y_var",
  title = "Top ROIs by |Age Beta|",
  plot_spec = list(theme = "minimal")
)
ggplot2::ggsave(file.path(out_root, "figures", "top_rois.png"), plot = top_rois_plot, width = 6, height = 5, bg = "white")

scatter_plot <- plot_scatter_assoc(
  df_combat,
  roi = scatter_roi,
  x = age_col,
  covariates = c(sex_col),
  color = sex_col,
  value = scatter_value,
  plot_spec = list(theme = "minimal", point_alpha = 0.55,
                   title = sprintf("%s vs %s (%s)", scatter_roi, age_col, scatter_value))
)
ggplot2::ggsave(file.path(out_root, "figures", "scatter_roi.png"), plot = scatter_plot, width = 6, height = 4, bg = "white")

# Optional brain map of Age beta (thresholded on FDR p <= 0.05)
brain_plot <- plot_brain_map_results(
  res_tbl,
  roi_col = "y_var",
  value_col = "estimate",
  atlas = "dk",
  roi_drop_suffix = roi_drop_suffix,
  plot_spec = list(
    layout = "dispersed",
    p_col = "p_adj",
    p_max = 0.05,
    legend_title = "Beta (Age)",
    title = "Age Association (FDR <= 0.05)",
    limit = NULL,
    ggseg_theme = "brain2"
  )
)
ggplot2::ggsave(file.path(out_root, "figures", "brain_map_age_beta.png"), plot = brain_plot, width = 8, height = 4, bg = "white")

t1 <- Sys.time()
timings$plots <- list(start = as.character(t0), end = as.character(t1),
                      seconds = as.numeric(difftime(t1, t0, units = "secs")))

## Step 5: Manifest + Logs -----------------------------------------------------
cat("\n## Step 5: Manifest + Logs\n")
t0 <- Sys.time()

cfg <- list(
  io = list(input = input_file, out_root = out_root),
  columns = list(site = site_col),
  roi = list(regex = roi_regex, drop_suffix = roi_drop_suffix),
  resolved = list(roi_cols = roi_cols),
  manifest = list(enabled = TRUE, path = file.path("logs", "run_manifest.json"))
)

ctx <- list(
  cfg = cfg,
  input = list(path = input_file, data = df0),
  data = df_combat,
  specs = list(prep = prep_spec_obj, combat = combat_spec_obj, analysis = analysis_spec_obj),
  results = list(assoc_scan = res_assoc),
  log = c(prep_res$log, combat_fit_obj$log, combat_apply_res$log, res_assoc$log),
  timings = timings,
  meta = list(
    batch_levels = levels(factor(df_combat[[site_col]]))
  )
)

manifest <- build_manifest(ctx)
write_manifest(manifest, file.path(out_root, "logs", "run_manifest.json"), overwrite = TRUE)
jsonlite::write_json(ctx$log, file.path(out_root, "logs", "run_log.json"), auto_unbox = TRUE, pretty = TRUE)

t1 <- Sys.time()
timings$manifest <- list(start = as.character(t0), end = as.character(t1),
                         seconds = as.numeric(difftime(t1, t0, units = "secs")))

## Summary ---------------------------------------------------------------------
cat("\n=== Association Demo Summary (Human) ===\n")
cat(sprintf("Input: %s\n", input_file))
cat(sprintf("Output root: %s\n", out_root))
cat(sprintf("Subjects: %d\n", nrow(df_combat)))
cat(sprintf("ROIs: %d\n", length(roi_cols)))
cat(sprintf("Min FDR p: %.4g\n", min(res_tbl$p_adj, na.rm = TRUE)))
top5 <- res_tbl[order(abs(res_tbl$estimate), decreasing = TRUE), c("y_var", "estimate", "p_adj")]
top5 <- head(top5, 5)
cat("Top 5 ROIs by |Age beta|:\n")
print(top5)
cat(sprintf("Scatter ROI: %s (%s)\n", scatter_roi, scatter_value))

t1_all <- Sys.time()
cat(sprintf("Total seconds: %.2f\n", as.numeric(difftime(t1_all, t0_all, units = "secs"))))
