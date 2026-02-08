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

# Human demo (multi-group):
# Step-by-step, explicit parameters, readable console output.
#
# NOTE: data/sample.csv is binary (Control/Case). For this multi-group demo, we
# deterministically split the Case group into two "subtypes" using the median
# Age among cases. In real usage, replace this with your actual subtype column.

# ------------------------------ Parameters -----------------------------------
input_file <- file.path("data", "sample.csv")
out_root <- file.path("outputs", "multiple_group_human")

# Columns (hard-coded for the sample dataset)
site_col <- "Site"
group_col_raw <- "Group"   # binary in sample.csv
group_col <- "Group3"      # derived 3-level group for this demo
age_col <- "Age"
sex_col <- "Sex"
icv_col <- "ICV"

# ROI detection
roi_regex <- "^(L_|R_)"
roi_drop_suffix <- "_thickavg"

# Group contrast direction: target - ref
ref_level <- "Control"
case_level <- "Case"
subtype_a <- "Case_A"
subtype_b <- "Case_B"

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

required <- c(site_col, group_col_raw, age_col, sex_col, icv_col)
missing <- setdiff(required, names(df0))
if (length(missing) > 0) stop(sprintf("Missing required columns: %s", paste(missing, collapse = ", ")), call. = FALSE)

roi_cols <- grep(roi_regex, names(df0), value = TRUE)
if (length(roi_cols) == 0) stop(sprintf("No ROI columns found with roi_regex='%s'.", roi_regex), call. = FALSE)

## Step 1: Prep ----------------------------------------------------------------
cat("\n## Step 1: Prep\n")
t0 <- Sys.time()

df0[[site_col]] <- factor(as.character(df0[[site_col]]))
df0[[sex_col]] <- factor(as.character(df0[[sex_col]]))

# Normalize case-control encodings (0/1 -> Control/Case)
df0[[group_col_raw]] <- {
  x <- df0[[group_col_raw]]
  if (all(unique(x) %in% c(0, 1, "0", "1"))) {
    factor(as.character(x), levels = c("0", "1"), labels = c(ref_level, case_level))
  } else {
    factor(as.character(x))
  }
}

if (!all(c(ref_level, case_level) %in% levels(df0[[group_col_raw]]))) {
  stop(sprintf("Expected %s to contain levels '%s' and '%s'. Found: %s",
               group_col_raw, ref_level, case_level,
               paste(levels(df0[[group_col_raw]]), collapse = ", ")),
       call. = FALSE)
}

# Deterministically split the Case group into 2 subtypes using median Age among cases.
is_case <- as.character(df0[[group_col_raw]]) == case_level
med_age <- stats::median(df0[[age_col]][is_case], na.rm = TRUE)
if (!is.finite(med_age)) stop("Cannot derive subtypes: median case Age is not finite.", call. = FALSE)

df0[[group_col]] <- ref_level
df0[[group_col]][is_case & df0[[age_col]] <= med_age] <- subtype_a
df0[[group_col]][is_case & df0[[age_col]] > med_age] <- subtype_b
df0[[group_col]] <- factor(df0[[group_col]], levels = c(ref_level, subtype_a, subtype_b))

sex_ref <- if ("0" %in% levels(df0[[sex_col]])) "0" else levels(df0[[sex_col]])[1]
factor_levels <- list()
factor_levels[[sex_col]] <- sex_ref
factor_levels[[group_col]] <- ref_level

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

## Step 3: Multi-Group Comparison ---------------------------------------------
cat("\n## Step 3: Compare Groups (Multi-Group)\n")
t0 <- Sys.time()

analysis_spec_obj <- analysis_spec(
  group_col = group_col,
  covariates = c(age_col, sex_col, icv_col),
  model_engine = "lm",
  diagnostics = "light",
  effect_size = "hedges_g_pooled",
  effect_size_ci = "wald",
  # For multi-group, compute effect sizes per contrast vs reference (target - ref).
  effect_size_scope = "pairwise",
  effect_size_from = "raw_data",
  p_adjust = list(method = "fdr", scope = "within_call",
                  family_id = "ct_multi_vs_control",
                  family_desc = "CT: subtype vs control"),
  ref_levels = setNames(list(ref_level), group_col),
  contrast = list(var = group_col, ref = ref_level, target = NULL, direction = "target_minus_ref")
)

res_comp <- compare_groups(df_combat, vars = roi_cols, spec = analysis_spec_obj, return = "result")
res_tbl <- res_comp$data

export_results_tables(
  res_tbl,
  out_dir = file.path(out_root, "tables"),
  base_name = "ct_multi_vs_control",
  formats = c("csv"),
  style = "publication",
  overwrite = TRUE,
  roi_drop_suffix = roi_drop_suffix
)
if (!is.null(res_comp$summaries)) {
  write.csv(res_comp$summaries, file.path(out_root, "tables", "ct_group_summary.csv"), row.names = FALSE)
}
if (!is.null(res_comp$diagnostics)) {
  write.csv(res_comp$diagnostics, file.path(out_root, "tables", "ct_diagnostics.csv"), row.names = FALSE)
}

# Also produce fully pairwise comparisons (Control vs Case_A, Control vs Case_B, Case_A vs Case_B)
analysis_spec_pw <- analysis_spec(
  group_col = group_col,
  covariates = c(age_col, sex_col, icv_col),
  model_engine = "lm",
  diagnostics = "none",
  effect_size = "hedges_g_pooled",
  effect_size_ci = "wald",
  effect_size_scope = "pairwise",
  effect_size_from = "raw_data",
  p_adjust = list(method = "fdr", scope = "within_call",
                  family_id = "ct_pairwise",
                  family_desc = "CT: pairwise subgroup comparisons"),
  ref_levels = list(),
  contrast = list(var = group_col, ref = NULL, target = NULL, direction = "target_minus_ref")
)

res_pw <- pairwise_groups(df_combat, vars = roi_cols, spec = analysis_spec_pw, return = "result")
export_results_tables(
  res_pw$data,
  out_dir = file.path(out_root, "tables"),
  base_name = "ct_pairwise_groups",
  formats = c("csv"),
  style = "publication",
  overwrite = TRUE,
  roi_drop_suffix = roi_drop_suffix
)

t1 <- Sys.time()
timings$compare_groups <- list(start = as.character(t0), end = as.character(t1),
                               seconds = as.numeric(difftime(t1, t0, units = "secs")))

## Step 4: Plots ---------------------------------------------------------------
cat("\n## Step 4: Plots\n")
t0 <- Sys.time()

for (tg in c(subtype_a, subtype_b)) {
  sub_tbl <- res_tbl[res_tbl$target_level == tg & res_tbl$ref_level == ref_level, , drop = FALSE]
  if (nrow(sub_tbl) == 0) next
  p <- plot_brain_map_results(
    sub_tbl,
    roi_col = "var",
    value_col = "es_value",
    atlas = "dk",
    roi_drop_suffix = roi_drop_suffix,
    plot_spec = list(layout = "dispersed", limit = NULL, legend_title = "Hedges g",
                     title = sprintf("%s - %s (ComBat-harmonized)", tg, ref_level),
                     ggseg_theme = "brain2")
  )
  ggplot2::ggsave(file.path(out_root, "figures", sprintf("brain_map_%s_vs_%s.png", tg, ref_level)),
                  plot = p, width = 8, height = 4, bg = "white")
}

qc_plot <- plot_qc_pvalue_hist(
  res_tbl,
  p_col = "p_adj",
  title = "QC: FDR p-value distribution (multi vs control)",
  plot_spec = list(theme = "minimal")
)
ggplot2::ggsave(file.path(out_root, "figures", "qc_pvalue_hist.png"), plot = qc_plot, width = 6, height = 4, bg = "white")

# Quadrant plot: subtype A vs Control (x) vs subtype B vs Control (y)
sub_a_tbl <- res_tbl[res_tbl$target_level == subtype_a & res_tbl$ref_level == ref_level, c("var", "es_value")]
sub_b_tbl <- res_tbl[res_tbl$target_level == subtype_b & res_tbl$ref_level == ref_level, c("var", "es_value")]
colnames(sub_a_tbl)[2] <- "es_a"
colnames(sub_b_tbl)[2] <- "es_b"
quad_tbl <- merge(sub_a_tbl, sub_b_tbl, by = "var", all = FALSE)
axis_max <- max(abs(c(quad_tbl$es_a, quad_tbl$es_b)), na.rm = TRUE)
if (!is.finite(axis_max) || axis_max == 0) axis_max <- 1
quad_plot <- plot_quadrant(
  quad_tbl,
  x_col = "es_a",
  y_col = "es_b",
  label_col = "var",
  axis_max = axis_max,
  add_label = FALSE,
  plot_spec = list(theme = "minimal",
                   title = sprintf("Quadrant: %s vs %s (x) and %s vs %s (y)",
                                   subtype_a, ref_level, subtype_b, ref_level))
)
ggplot2::ggsave(file.path(out_root, "figures", "quadrant_subtypes_vs_control.png"),
                plot = quad_plot, width = 6, height = 5, bg = "white")

# Raincloud plot with pairwise significance markers for a selected ROI
rain_roi <- "L_superiorfrontal_thickavg"
rain_plot <- plot_raincloud_roi(
  df_combat,
  roi = rain_roi,
  group_col = group_col,
  covariates = c(age_col, sex_col, icv_col),
  value = "marginal",
  pairwise = list(method = "lm", p_adjust = "holm", label = "p_adj", p_col = "p_adj",
                  show_ns = TRUE),
  plot_spec = list(theme = "minimal", legend_position = "none",
                   title = sprintf("Raincloud: %s (marginal)", rain_roi))
)
ggplot2::ggsave(file.path(out_root, "figures", "raincloud_roi.png"),
                plot = rain_plot, width = 7, height = 4, bg = "white")

t1 <- Sys.time()
timings$plots <- list(start = as.character(t0), end = as.character(t1),
                      seconds = as.numeric(difftime(t1, t0, units = "secs")))

## Step 5: Manifest + Logs -----------------------------------------------------
cat("\n## Step 5: Manifest + Logs\n")
t0 <- Sys.time()

# Minimal cfg (only used for manifest builder)
cfg <- list(
  io = list(input = input_file, out_root = out_root),
  columns = list(site = site_col),
  roi = list(regex = roi_regex, drop_suffix = roi_drop_suffix),
  resolved = list(group_col = group_col, roi_cols = roi_cols),
  manifest = list(enabled = TRUE, path = file.path("logs", "run_manifest.json"))
)

ctx <- list(
  cfg = cfg,
  input = list(path = input_file, data = df0),
  data = df_combat,
  specs = list(prep = prep_spec_obj, combat = combat_spec_obj, analysis = analysis_spec_obj, pairwise = analysis_spec_pw),
  results = list(compare_groups = res_comp, pairwise_groups = res_pw),
  log = c(prep_res$log, combat_fit_obj$log, combat_apply_res$log, res_comp$log, res_pw$log),
  timings = timings,
  meta = list(
    group_counts = as.list(table(df_combat[[group_col]])),
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
cat("\n=== Multi-Group Demo Summary (Human) ===\n")
cat(sprintf("Input: %s\n", input_file))
cat(sprintf("Output root: %s\n", out_root))
cat(sprintf("Subjects: %d\n", nrow(df_combat)))
cat(sprintf("ROIs: %d\n", length(roi_cols)))
cat(sprintf("Group column: %s (derived from %s)\n", group_col, group_col_raw))
cat("Group counts:\n")
print(table(df_combat[[group_col]]))

cat(sprintf("Min FDR p: %.4g\n", min(res_tbl$p_adj, na.rm = TRUE)))

for (tg in c(subtype_a, subtype_b)) {
  sub_tbl <- res_tbl[res_tbl$target_level == tg & res_tbl$ref_level == ref_level, , drop = FALSE]
  if (nrow(sub_tbl) == 0) next
  top5 <- sub_tbl[order(abs(sub_tbl$es_value), decreasing = TRUE), c("var", "es_value", "p_adj")]
  top5 <- head(top5, 5)
  cat(sprintf("\nTop 5 ROIs by |effect size| (%s - %s):\n", tg, ref_level))
  print(top5)
}

t1_all <- Sys.time()
cat(sprintf("\nTotal seconds: %.2f\n", as.numeric(difftime(t1_all, t0_all, units = "secs"))))
