# roiflow API: Usage + Arguments

This reference is generated from `man/*.Rd` and shows the exact function signatures and documented parameters.

## `analysis_spec`

Create Analysis Specification

**Usage**

```r
analysis_spec(
  group_col = "groups",
  covariates = c("Age", "Sex", "ICV"),
  model_engine = c("auto", "lm", "glm", "lmer", "glmer"),
  family = NULL,
  random_effects = NULL,
  site_col = NULL,
  site_effect = c("none", "fixed", "random"),
  diagnostics = c("none", "light", "full"),
  effect_size = c("hedges_g_pooled", "cohen_d_pooled", "glass_delta",
    "standardized_beta", "none"),
  effect_size_ci = c("wald", "boot", "none"),
  effect_size_from = c("raw_data", "model"),
  effect_size_scope = c("two_group_only", "pairwise", "none"),
  p_adjust = list(method = "fdr", scope = "within_call", family_id = NA_character_,
    family_desc = NA_character_),
  conf_int = FALSE,
  conf_level = 0.95,
  ref_levels = list(),
  contrast = list(var = NULL, ref = NULL, target = NULL, direction = "target_minus_ref")
)
```

**Arguments**

- `group_col`: Character. Grouping column for group comparisons.
- `covariates`: Character vector. Covariate column names.
- `model_engine`: Character. One of "auto","lm","glm","lmer","glmer".
- `family`: Optional family for GLM/GLMER. Either a family object or a string such as "gaussian" or "binomial".
- `random_effects`: Character. Random effects term, e.g. "(1|Site)".
- `site_col`: Optional site column name.
- `site_effect`: Character. How to use site_col: "none", "fixed", or "random".
- `diagnostics`: Character. One of "none","light","full".
- `effect_size`: Character. Effect size type. Default is Hedges' g for two groups.
- `effect_size_ci`: Character. Effect size CI method: "wald" (default), "boot", or "none".
- `effect_size_from`: Character. "raw_data" or "model".
- `effect_size_scope`: Character. "two_group_only", "pairwise", or "none".
- `p_adjust`: List. Multiple-testing specification with fields: method, scope, family_id, family_desc.
- `conf_int`: Logical. Whether to compute confidence intervals for model terms.
- `conf_level`: Numeric. Confidence level.
- `ref_levels`: Named list. Reference levels for categorical variables (e.g., list(Group = "Control")).
- `contrast`: List. Contrast specification with fields var, ref, target, direction. Direction must be "target_minus_ref".

## `assoc_scan`

Association Scan Across Outcomes

**Usage**

```r
assoc_scan(
  df,
  y_vars,
  x,
  spec = analysis_spec(),
  by_group = FALSE,
  interaction = FALSE,
  return = c("result", "data")
)
```

**Arguments**

- `df`: Data frame.
- `y_vars`: Character vector of outcome variables.
- `x`: Character. Predictor variable.
- `spec`: Analysis specification.
- `by_group`: Logical. If TRUE, run separately within each group.
- `interaction`: Logical. If TRUE, include x * group_col.
- `return`: Character. "result" or "data".

## `build_manifest`

Build a Run Manifest From Pipeline Context

**Usage**

```r
build_manifest(ctx)
```

**Arguments**

- `ctx`: Pipeline context returned by run_pipeline.

## `combat_apply`

Apply a ComBat Harmonization Model

**Usage**

```r
combat_apply(tmp, dat, batch, mod = NULL, verbose = TRUE)
```

**Arguments**

- `tmp`: Fitted model list returned by combat_fit.
- `dat`: Numeric matrix or data frame. Either rows or columns must match the length of batch. If the fitted model stored a transpose flag, the same orientation is enforced.
- `batch`: Factor indicating batch/site per subject. Levels must match the levels used during fit.
- `mod`: Optional numeric design matrix of covariates. Must match the fit-time design (column count and names); missing fit-time columns are zero-padded, but new columns are rejected.
- `verbose`: Logical. Print progress messages.

## `combat_apply_df`

Apply ComBat to a Data Frame

**Usage**

```r
combat_apply_df(fit, df, verbose = interactive())
```

**Arguments**

- `fit`: Object returned by combat_fit_df.
- `df`: Data frame to harmonize.
- `verbose`: Logical. Print progress messages.

## `combat_fit`

Fit a ComBat Harmonization Model

**Usage**

```r
combat_fit(
  dat,
  batch,
  mod = NULL,
  eb = TRUE,
  verbose = TRUE,
  single_batch = c("error", "return")
)
```

**Arguments**

- `dat`: Numeric matrix or data frame. Either rows or columns must match the length of batch. If subjects are rows, the data are transposed internally. Features are assumed to be rows after any transpose.
- `batch`: Factor indicating batch/site per subject. Must be a factor; unused levels are dropped with a warning.
- `mod`: Optional numeric design matrix of covariates (no intercept required). Must have the same number of rows as length(batch) when provided.
- `eb`: Logical. If TRUE, use empirical Bayes shrinkage (default).
- `verbose`: Logical. Print progress messages.
- `single_batch`: Character. What to do if there is only one batch level: "error" (default) stops; "return" returns a no-op model that will pass data through unchanged in combat_apply.

## `combat_fit_df`

Fit ComBat on a Data Frame

**Usage**

```r
combat_fit_df(df, spec = combat_spec(), verbose = interactive())
```

**Arguments**

- `df`: Data frame with ROI columns and batch column.
- `spec`: Specification list from combat_spec.
- `verbose`: Logical. Print progress messages.

## `combat_one_liner`

One-Liner for Fit + Apply

**Usage**

```r
combat_one_liner(
  df,
  spec = combat_spec(),
  return = c("data", "result"),
  verbose = interactive()
)
```

**Arguments**

- `df`: Data frame with ROI columns and batch column.
- `spec`: Specification list from combat_spec.
- `return`: Character. "data" returns only the harmonized data frame; "result" returns the full combat_result object.
- `verbose`: Logical. Print progress messages.

## `combat_spec`

Define ComBat Wrapper Specifications

**Usage**

```r
combat_spec(
  roi_cols = NULL,
  roi_regex = "^(L_|R_)",
  batch_col = "Site",
  mod_formula = NULL,
  eb = TRUE,
  missing_fit = c("impute", "error"),
  single_batch = c("error", "noop"),
  keep_inputs = TRUE
)
```

**Arguments**

- `roi_cols`: Character vector. Explicit ROI column names. If NULL, columns are inferred via roi_regex.
- `roi_regex`: Character. Regex to identify ROI columns when roi_cols is NULL. Default matches L_/R_ prefixed columns.
- `batch_col`: Character. Batch/site column name.
- `mod_formula`: Model formula for covariates (e.g., ~ Age + Sex). Use NULL for no covariates.
- `eb`: Logical. Use empirical Bayes shrinkage in ComBat.
- `missing_fit`: Character. How to handle missing values during fit: "impute" (default) or "error" to stop if any NA is present in ROI data.
- `single_batch`: Character. If only one batch level is present: "error" (default) or "noop" to return a pass-through model.
- `keep_inputs`: Logical. If TRUE, return all columns on apply; otherwise return only ROI columns.

## `compare_groups`

Compare Groups Across Variables

**Usage**

```r
compare_groups(df, vars, spec = analysis_spec(), return = c("result", "data"))
```

**Arguments**

- `df`: Data frame with outcomes, group column, and covariates.
- `vars`: Character vector of outcome variables.
- `spec`: Analysis specification from analysis_spec.
- `return`: Character. "result" (default) returns a structured result object; "data" returns only the results table.

## `corr_scan`

Correlation Scan

**Usage**

```r
corr_scan(
  df,
  vars,
  x,
  method = c("pearson", "spearman"),
  covariates = NULL,
  return = c("result", "data")
)
```

**Arguments**

- `df`: Data frame.
- `vars`: Character vector of variables to correlate with x.
- `x`: Character. Target variable.
- `method`: Character. "pearson" or "spearman".
- `covariates`: Optional covariates for partial correlation.
- `return`: Character. "result" or "data".

## `export_results_tables`

Export Raw + Formatted Results Tables

**Usage**

```r
export_results_tables(
  res_tbl,
  out_dir,
  base_name = "results",
  formats = c("csv"),
  style = c("publication", "minimal"),
  overwrite = FALSE,
  ...
)
```

**Arguments**

- `res_tbl`: Data frame. Raw results table.
- `out_dir`: Character. Output directory for tables.
- `base_name`: Character. File prefix (without extension).
- `formats`: Character vector. Any of "csv" or "tsv".
- `style`: Character. Formatting style for the formatted table.
- `overwrite`: Logical. If FALSE (default), error if outputs exist.
- `...`: Passed to format_results_table.

## `format_results_table`

Format a Results Table

**Usage**

```r
format_results_table(
  res_tbl,
  style = c("minimal", "publication"),
  roi_drop_suffix = "_thickavg",
  roi_drop_prefix = "^SV_",
  digits = 3,
  p_digits = 3
)
```

**Arguments**

- `res_tbl`: Data frame. Standard results table.
- `style`: Character. One of "minimal" or "publication".
- `roi_drop_suffix`: Optional character. If provided, remove this suffix from ROI names (e.g., "_thickavg").
- `roi_drop_prefix`: Optional character. Regex prefix to remove from ROI names (e.g., "^SV_").
- `digits`: Integer. Rounding digits for effect/estimate columns.
- `p_digits`: Integer. Significant digits for p-values.

## `make_report`

Save Results, Plots, and Manifest

**Usage**

```r
make_report(result, out_dir, report_spec = list())
```

**Arguments**

- `result`: A profile_result object.
- `out_dir`: Output directory.
- `report_spec`: List with options: table_format (csv/tsv), plot_format (png/pdf), write_manifest (TRUE/FALSE).

## `pairwise_groups`

Pairwise Group Comparisons

**Usage**

```r
pairwise_groups(df, vars, spec = analysis_spec(), return = c("result", "data"))
```

**Arguments**

- `df`: Data frame with outcomes and group column.
- `vars`: Character vector of outcome variables.
- `spec`: Analysis specification from analysis_spec.
- `return`: Character. "result" (default) or "data".

## `pairwise_tests_roi`

Pairwise Tests for a Single ROI Across Groups

**Usage**

```r
pairwise_tests_roi(
  data,
  roi,
  group_col,
  covariates = NULL,
  value = c("raw", "marginal"),
  method = c("t_test", "wilcox", "lm"),
  p_adjust = c("holm", "fdr", "bonferroni", "none"),
  pairs = NULL
)
```

**Arguments**

- `data`: Data frame.
- `roi`: Character. ROI column name.
- `group_col`: Character. Grouping column (>= 2 levels).
- `covariates`: Optional character vector of covariate column names.
- `value`: Character. "raw" or "marginal". Only used for t_test/wilcox (marginal values are derived from ROI ~ group + covariates by removing covariate contributions while keeping group effects).
- `method`: Character. One of "t_test", "wilcox", "lm".
- `p_adjust`: Character. Passed to p.adjust. Use "none" to skip adjustment.
- `pairs`: Optional list of length-2 character vectors specifying which pairs to test. Default: all pairs.

## `plot_auc_hist`

AUC Histogram Plot

**Usage**

```r
plot_auc_hist(
  auc_null,
  auc_true = NULL,
  nbins = 30,
  xlim = NULL,
  ylim = NULL,
  line_color = "steelblue1",
  plot_spec = list()
)
```

**Arguments**

- `auc_null`: Numeric vector. Null AUC distribution.
- `auc_true`: Optional numeric. Observed/true AUC.
- `nbins`: Integer. Histogram bin count.
- `xlim`: Optional numeric vector of length 2. X limits.
- `ylim`: Optional numeric vector of length 2. Y limits.
- `line_color`: Character. Color for the observed AUC line.
- `plot_spec`: List. Optional plot customization (theme, fill/alpha).

## `plot_bar`

Bar Plot of Group Profiles

**Usage**

```r
plot_bar(
  summary_tbl,
  se_col = "se",
  group_col = "groups",
  var_col = "var",
  value_col = "mean",
  plot_spec = list()
)
```

**Arguments**

- `summary_tbl`: Data frame with summary statistics.
- `se_col`: Character. Column name for standard error.
- `group_col`: Character. Column name for group labels.
- `var_col`: Character. Column name for variable/ROI names.
- `value_col`: Character. Column name for the mean (or plotted value).
- `plot_spec`: List. Optional plot customization (theme, sizes, colors).

## `plot_brain_map`

Brain Map Plot (ggseg-based)

**Usage**

```r
plot_brain_map(tbl, atlas_spec, plot_spec = list())
```

**Arguments**

- `tbl`: Data frame with at least a label column (e.g. "lh_bankssts").
- `atlas_spec`: List with fields: atlas (dk/desterieux/aseg), label_col, value_col.
- `plot_spec`: List. Common fields: layout: "dispersed" (default) or "stacked". view: optional view filter (e.g. "lateral"). hemisphere: optional hemi filter (e.g. c("left","right")). p_col, p_max: optional thresholding (set fill to NA if p > p_max). limit, midpoint, legend_title, title.

## `plot_brain_map_results`

Brain Map from a Results Table (ROI Names -> ggseg Labels)

**Usage**

```r
plot_brain_map_results(
  res_tbl,
  roi_col = "var",
  value_col = "es_value",
  atlas = c("dk", "desterieux", "aseg"),
  roi_drop_suffix = "_thickavg",
  plot_spec = list()
)
```

**Arguments**

- `res_tbl`: Data frame (typically a results table).
- `roi_col`: Character. Column containing ROI names. Default: "var".
- `value_col`: Character. Column containing the value to plot.
- `atlas`: Character. One of "dk", "desterieux", "aseg".
- `roi_drop_suffix`: Optional suffix removed from ROI names before mapping.
- `plot_spec`: List. Passed through to plot_brain_map.

## `plot_dot`

Dot Plot of Group Profiles

**Usage**

```r
plot_dot(
  summary_tbl,
  se_col = "se",
  group_col = "groups",
  var_col = "var",
  value_col = "mean",
  plot_spec = list()
)
```

**Arguments**

- `summary_tbl`: Data frame with summary statistics (e.g. from compare_groups()$summaries).
- `se_col`: Character. Column name for standard error.
- `group_col`: Character. Column name for group labels.
- `var_col`: Character. Column name for variable/ROI names.
- `value_col`: Character. Column name for the mean (or plotted value).
- `plot_spec`: List. Optional plot customization (theme, sizes, colors).

## `plot_qc_pvalue_hist`

QC Plot: P-Value Histogram

**Usage**

```r
plot_qc_pvalue_hist(
  res_tbl,
  p_col = "p_adj",
  nbins = 30,
  title = "QC: p-value distribution",
  plot_spec = list()
)
```

**Arguments**

- `res_tbl`: Data frame (results table).
- `p_col`: Character. Column name containing p-values (e.g. "p_adj").
- `nbins`: Integer. Histogram bin count.
- `title`: Character. Plot title.
- `plot_spec`: List. Optional plot customization (theme, fill/alpha).

## `plot_quadrant`

Quadrant Plot from Stat Table

**Usage**

```r
plot_quadrant(
  stat_tbl,
  x_col,
  y_col,
  label_col = NULL,
  axis_max = 1,
  add_label = FALSE,
  plot_spec = list()
)
```

**Arguments**

- `stat_tbl`: Data frame (typically a results table).
- `x_col`: Character. Column name for x-axis value.
- `y_col`: Character. Column name for y-axis value.
- `label_col`: Optional character. Column name used for point labels.
- `axis_max`: Numeric. Symmetric axis limit (±axis_max).
- `add_label`: Logical. Add repelled labels.
- `plot_spec`: List. Optional plot customization (theme, point size/alpha).

## `plot_radar`

Radar Plot for Group Profiles

**Usage**

```r
plot_radar(
  summary_tbl,
  var_col = "var",
  value_col = "mean",
  group_col = "groups",
  vars = NULL,
  max_value = 2,
  min_value = -2,
  plot_spec = list()
)
```

**Arguments**

- `summary_tbl`: Data frame with summary statistics.
- `var_col`: Character. Column name for variable/ROI names.
- `value_col`: Character. Column name for the mean (or plotted value).
- `group_col`: Character. Column name for group labels.
- `vars`: Optional character vector to subset variables.
- `max_value`: Numeric. Max axis value used by radar chart.
- `min_value`: Numeric. Min axis value used by radar chart.
- `plot_spec`: List. Optional plot customization (theme, title, base_size).

## `plot_raincloud_roi`

Raincloud Plot for a Selected ROI (Half-Violin + Points)

**Usage**

```r
plot_raincloud_roi(
  data,
  roi,
  group_col,
  covariates = NULL,
  value = c("raw", "marginal"),
  orientation = c("vertical", "horizontal"),
  style = c("half_violin"),
  pairwise = NULL,
  plot_spec = list()
)
```

**Arguments**

- `data`: Data frame.
- `roi`: Character. ROI column name.
- `group_col`: Character. Grouping column (>= 2 levels).
- `covariates`: Optional character vector of covariate column names.
- `value`: Character. "raw" or "marginal".
- `orientation`: Character. "vertical" (default) or "horizontal".
- `style`: Character. Currently only "half_violin" is supported.
- `pairwise`: Optional list enabling pairwise comparisons + markers. Fields: method: "t_test", "wilcox", or "lm". p_adjust: "holm", "fdr", "bonferroni", "none". pairs: optional list of pairs (each length-2 character vector). label: "stars" (default) or "p_adj". p_col: "p_adj" (default) or "p_value". show_ns: show non-significant labels (default: FALSE).
- `plot_spec`: List. Customization options (theme, sizes, colors, etc.).

## `plot_scatter_assoc`

Scatter Plot for Association (Selected ROI)

**Usage**

```r
plot_scatter_assoc(
  data,
  roi,
  x,
  covariates = NULL,
  color = NULL,
  value = c("raw", "marginal"),
  plot_spec = list()
)
```

**Arguments**

- `data`: Data frame.
- `roi`: Character. ROI column name.
- `x`: Character. Predictor variable (x-axis).
- `covariates`: Optional character vector of covariates to adjust for.
- `color`: Optional character. Column used for point colors (e.g. Sex).
- `value`: Character. "raw" or "marginal".
- `plot_spec`: List. Customization options (theme, sizes, colors).

## `plot_top_rois`

Top ROIs Bar Plot

**Usage**

```r
plot_top_rois(
  res_tbl,
  metric = "es_value",
  n = 20,
  label_col = "var",
  title = "Top ROIs",
  plot_spec = list()
)
```

**Arguments**

- `res_tbl`: Data frame (results table).
- `metric`: Character. Column used for ranking/plotting (absolute value).
- `n`: Integer. Number of rows to plot.
- `label_col`: Character. Column used as labels.
- `title`: Character. Plot title.
- `plot_spec`: List. Optional plot customization (theme, gradient colors).

## `plot_violin`

Violin Plot (Raw Data, No Re-fitting)

**Usage**

```r
plot_violin(data, x, y, group_col = NULL, plot_spec = list())
```

**Arguments**

- `data`: Data frame.
- `x`: Character. X-axis variable (typically a grouping column).
- `y`: Character. Y-axis variable.
- `group_col`: Optional character. Column used for fill/group aesthetics.
- `plot_spec`: List. Optional plot customization (theme, sizes, colors).

## `prep`

Preprocess ROI Data

**Usage**

```r
prep(
  df,
  spec = prep_spec(),
  return = c("data", "result"),
  verbose = interactive(),
  strict = FALSE
)
```

**Arguments**

- `df`: Data frame containing ROI measurements and covariates.
- `spec`: Preprocessing specification created by prep_spec. Default uses standard settings for L_/R_ prefixed ROI columns.
- `return`: Character. What to return: "data": Return only the cleaned data frame (default) "result": Return full result object with data, QC metrics, and logs
- `verbose`: Logical. Print progress summary. Default: TRUE in interactive sessions, FALSE otherwise.
- `strict`: Logical. If TRUE, stop on warnings (e.g., missing reference levels, no ROI columns). Default: FALSE.

## `prep_spec`

Create Preprocessing Specification

**Usage**

```r
prep_spec(
  roi_regex = "^(L_|R_)",
  numeric_regex = "^(L_|R_)|Age$|ICV$|mod_PDS$",
  zero_to_na_regex = "^(L_|R_)|^ICV$",
  factor_levels = list(Sex = "0", Group = "0"),
  site_col = "Site",
  site_ref = c("largest", "keep"),
  icv_col = "ICV",
  icv_scale = 1e+06,
  na_action = c("drop_required", "drop_all", "keep"),
  required_cols = NULL,
  outlier_action = c("clip", "flag", "none"),
  outlier_method = c("iqr", "sd"),
  outlier_k = 3
)
```

**Arguments**

- `roi_regex`: Character. Regular expression to identify ROI columns. Default: "^(L_|R_)" matches columns starting with L_ or R_.
- `numeric_regex`: Character. Regular expression to identify columns that should be converted to numeric. Default matches ROI columns, Age, ICV, and mod_PDS.
- `zero_to_na_regex`: Character. Regular expression to identify columns where zeros should be converted to NA. Default matches ROI columns and ICV.
- `factor_levels`: Named list. Specifies factor variables and their reference levels. Default: list(Sex = "0", Group = "0").
- `site_col`: Character. Name of the site/scanner column. Default: "Site".
- `site_ref`: Character or vector. How to set the site reference level: "largest": Use the site with most observations (default) "keep": Keep existing factor order A specific site name: Use that site as reference
- `icv_col`: Character. Name of the intracranial volume column. Default: "ICV".
- `icv_scale`: Numeric. Factor to divide ICV by. Default: 1e6 (converts to millions).
- `na_action`: Character. How to handle missing data: "drop_required": Drop rows with NA in required columns (default) "drop_all": Drop rows with any NA "keep": Keep all rows
- `required_cols`: Character vector. Columns that must be complete for na_action = "drop_required". If NULL (default), inferred as key covariates (Site, Sex, Group, Age, ICV if present).
- `outlier_action`: Character. How to handle outliers: "clip": Clip outliers to fence boundaries (default) "flag": Detect but don't modify "none": Skip outlier detection
- `outlier_method`: Character. Method for outlier detection: "iqr": Interquartile range method (default) "sd": Standard deviation method
- `outlier_k`: Numeric. Multiplier for outlier fences. Default: 3. For IQR: fences at Q1 - k*IQR and Q3 + k*IQR. For SD: fences at mean ± k*SD.

## `print.combat_result`

Print Method for combat_result

**Usage**

```r
list(list("print"), list("combat_result"))(x, ...)
```

**Arguments**

- `x`: Object of class combat_result.
- `...`: Unused.

## `print.roiflow_prep_result`

Print Method for Preprocessing Results

**Usage**

```r
list(list("print"), list("roiflow_prep_result"))(x, ...)
```

**Arguments**

- `x`: Object of class roiflow_prep_result.
- `...`: Additional arguments (currently unused).

## `read_cfg`

Read a Pipeline Config (JSON)

**Usage**

```r
read_cfg(path)
```

**Arguments**

- `path`: Character. Path to a JSON file.

## `resolve_cfg`

Resolve a Config Against a Data Frame

**Usage**

```r
resolve_cfg(cfg, df = NULL)
```

**Arguments**

- `cfg`: Config list (validated).
- `df`: Optional data frame. If provided, enables auto-detection.

## `roiflow`

Package-Oriented Alias for run_pipeline

**Usage**

```r
roiflow(cfg_or_path)
```

**Arguments**

- `cfg_or_path`: Config list or JSON path.

## `run_pipeline`

Run a Config-First Analysis Pipeline

**Usage**

```r
run_pipeline(cfg_or_path)
```

**Arguments**

- `cfg_or_path`: Either a config list or a path to a JSON config file.

## `validate_cfg`

Validate and Normalize a Pipeline Config

**Usage**

```r
validate_cfg(cfg)
```

**Arguments**

- `cfg`: A nested list.

## `write_manifest`

Write a Manifest to JSON

**Usage**

```r
write_manifest(manifest, path, overwrite = FALSE)
```

**Arguments**

- `manifest`: List from build_manifest.
- `path`: Output file path.
- `overwrite`: Logical. If FALSE, error if file exists.

