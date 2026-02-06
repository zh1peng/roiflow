# roiflow Job JSON Template (Config-First Pipeline)

This document explains the JSON config schema used by:

- `run_pipeline(cfg_or_path)`
- `roiflow(cfg_or_path)` (alias)

The template file lives at `inst/examples/job_template.json`.

## Run It

```r
devtools::load_all(".")
ctx <- run_pipeline("inst/examples/job_template.json")
```

All paths are interpreted relative to your current working directory (usually the project root when running `Rscript`).

## Top-Level Keys

- `project`: Optional metadata (`name`, `version`) used in the manifest.
- `io`: Input/output and saving behavior.
- `columns`: Column mapping (site, group, covariates).
- `roi`: ROI column detection.
- `pipeline`: Ordered list of steps.
- `manifest`: Whether/where to write `run_manifest.json`.

## `io`

- `io.input` (string, required): Path to CSV input.
- `io.out_root` (string, required): Output root folder. The pipeline writes into `data/`, `tables/`, `figures/`, `logs/`.
- `io.overwrite` (bool): If `false`, the pipeline errors if an output file already exists.
- `io.save_data` (bool): Write the final dataset to `out_root/data/`.
- `io.save_tables` (bool): Write results tables to `out_root/tables/`.
- `io.save_plots` (bool): Write plots to `out_root/figures/`.
- `io.plot_formats` (string array): Formats passed to `ggplot2::ggsave()`. Default is `["png"]`.

## `columns`

- `columns.site` (string): Site/batch column for ComBat and optional diagnostics/manifest info.
- `columns.group` (string): Group column for group comparisons.
  - Use `"auto"` to auto-detect from common names (e.g., `Group`, `Diagnosis`).
- `columns.age` / `columns.sex` / `columns.icv` (string or `null`): Convenience names for common covariates.
  - If your dataset uses different covariate names, set them here (or pass explicit covariate lists per step).

## `roi`

- `roi.regex` (string): Regex used to infer ROI columns when `roi.cols` is `null`.
- `roi.drop_suffix` (string): Optional suffix removed when formatting ROI labels (tables/plots).
- `roi.cols` (string array or `null`): Optional explicit ROI column list (bypasses regex detection).

## `pipeline`

An array of `{ "step": "<name>", "params": { ... } }`.

Supported step names:

- `prep`
- `combat`
- `compare_groups` (two-group only; requires `ref` and `target`)
- `pairwise_groups` (multi-group safe; runs all pairwise contrasts)
- `assoc_scan`
- `corr_scan`
- `plots`
- `export`

### Step: `prep`

Creates and runs `prep_spec()` + `prep()`.

Parameters (`pipeline[].params`):

- `group_ref` (string): Reference level for the resolved group column.
- `sex_ref` (string): Reference level for the `Sex` column (if present).
- `factor_levels` (object or `null`): Named mapping of factor column -> reference level. If provided, these are used; `group_ref` is still enforced if set.
- `na_action` (string): One of `drop_required`, `drop_all`, `keep`. Passed to `prep_spec()`.
- `outlier_action` (string): One of `none`, `flag`, `clip`. Passed to `prep_spec()`.
- `verbose` (bool): Passed to `prep()`.

### Step: `combat`

Fits and applies ComBat on ROI columns using `combat_fit_df()` + `combat_apply_df()`.

Parameters:

- `covariates` (string array): Covariates used in the ComBat design (becomes `~ Age + Sex + ICV`). Use `[]` for no covariates.
- `missing_fit` (string): `impute` (default) or `error`.
- `single_batch` (string): `error` (default) or `noop` (pass-through when only one batch exists).
- `verbose` (bool): Print ComBat progress messages.

### Step: `compare_groups` (two-group)

Runs `analysis_spec()` + `compare_groups()` on all ROIs.

Parameters:

- `covariates` (string array): Covariates to include in the model.
- `ref` (string, required): Reference level for the group factor.
- `target` (string, required): Target level for the group contrast. Direction is always `target - ref`.
- `model_engine` (string): One of `lm`, `glm`, `lmer`, `glmer`, `auto`.
- `diagnostics` (string): One of `none`, `light`, `full`.
- `effect_size` (string): One of `hedges_g_pooled`, `cohen_d_pooled`, `glass_delta`, `standardized_beta`, `none`.
- `effect_size_ci` (string): `wald` (default), `boot`, or `none`.
- `effect_size_from` (string): `raw_data` or `model`.
- `p_adjust_method` (string): `fdr` (default), or any method accepted by `p.adjust()`.
- `family_id` (string): Stored into results metadata as the multiple-testing family ID.
- `family_desc` (string): Stored into results metadata as the multiple-testing family description.

### Step: `pairwise_groups` (multi-group)

Runs `pairwise_groups()` for all group pairs using the same schema as `compare_groups()`.

Parameters:

- `covariates`, `model_engine`, `diagnostics`, `effect_size`, `p_adjust_method`, `family_id`, `family_desc` as above.
- `effect_size_ci` (string): `wald` (default), `boot`, or `none`.
- `ref_levels` (object): Optional reference levels for other categorical variables (rare).

### Step: `assoc_scan`

Association scan over ROIs: fits `ROI ~ x + covariates` (and optional group interaction).

Parameters:

- `x` (string, required): Predictor variable of interest.
- `covariates` (string array): Covariates to include.
- `model_engine`, `diagnostics`, `effect_size`, `p_adjust_method`, `family_id`, `family_desc` as above.
- `by_group` (bool): If `true`, run association separately within each group.
- `interaction` (bool): If `true`, fit `ROI ~ x * group + covariates`.

### Step: `corr_scan`

Correlation scan over ROIs.

Parameters:

- `x` (string, required): Predictor variable.
- `method` (string): `pearson` (default) or `spearman`.
- `covariates` (string array): Optional covariates for partial correlation via residualization.

### Step: `plots`

Builds ggplot objects from the results tables (no refitting).

Parameters:

- `enabled` (string array): Any of `brain_map`, `pvalue_hist`, `top_rois`, `raincloud_roi`, `scatter_assoc`.
- `brain_layout` (string): Passed to `plot_brain_map()` (e.g., `"dispersed"`).
- `atlas` (string): `"dk"` (default). Requires `ggseg` installed.
- `stat_col` (string): Column to plot on the brain map (default `es_value`).
- `p_type` (string): Column used for brain-map thresholding (default `p_adj`).
- `p_max` (number or `null`): If set, values with `p_type > p_max` are masked (`NA`) on the brain map.
- `p_col` (string): Column for p-value histogram (default `p_adj`).
- `title`, `legend_title`, `na_color` (strings): Plot appearance options (brain map).
- `ggseg_theme` (string): Brain-map theme (ggseg-native): `brain`, `brain2` (default), `darkbrain`, `custombrain`.
- `theme` (string): Optional global theme applied to non-brain plots (e.g., `minimal`, `classic`, `bw`).

Raincloud ROI (optional):

- `raincloud_roi.roi` (string, required if enabled): ROI column name.
- `raincloud_roi.value` (string): `raw` or `marginal` (covariate-adjusted).
- `raincloud_roi.covariates` (string array): Covariates to adjust out when `value="marginal"`.
- `raincloud_roi.group_col` (string or `null`): Group column (defaults to resolved `columns.group`).
- `raincloud_roi.orientation` (string): `vertical` (default) or `horizontal`.
- `raincloud_roi.style` (string): `half_violin` (default).
- `raincloud_roi.pairwise` (object or `null`): Pairwise comparisons + markers.
- `raincloud_roi.pairwise.method` (string): `t_test`, `wilcox`, or `lm` (default).
- `raincloud_roi.pairwise.p_adjust` (string): `holm` (default), `fdr`, `bonferroni`, `none`.
- `raincloud_roi.pairwise.pairs` (array or `null`): Explicit pairs to display (each item is a 2-string array).
- `raincloud_roi.pairwise.label` (string): `stars` (default) or `p_adj` (show adjusted p-values).
- `raincloud_roi.pairwise.p_col` (string): `p_adj` (default) or `p_value`.
- `raincloud_roi.pairwise.show_ns` (boolean): Show non-significant pairs (default: false).
- `raincloud_roi.plot_spec` (object): Passed to `plot_raincloud_roi()` (theme, sizes, colors, etc.).

Scatter association (optional):

- `scatter_assoc.roi` (string, required if enabled): ROI column name.
- `scatter_assoc.x` (string, required if enabled): Predictor variable (e.g., `Age`).
- `scatter_assoc.value` (string): `raw` or `marginal` (adjust for other covariates).
- `scatter_assoc.covariates` (string array): Covariates to adjust out when `value="marginal"`.
- `scatter_assoc.color` (string or `null`): Optional column used for point colors (e.g., `Sex`).
- `scatter_assoc.plot_spec` (object): Passed to `plot_scatter_assoc()` (theme, point alpha/size, etc.).

### Step: `export`

Writes datasets, tables, plots, logs, and manifest under `io.out_root/`.

Parameters:

- `base_name` (string): Base name for the main `compare_groups` results (default `ct_group_diff`).
- `formats` (string array): Any of `csv` or `tsv`.
- `style` (string): `publication` (default) or `minimal` formatting for results tables.

## `manifest`

- `manifest.enabled` (bool): If `true`, writes a derived manifest after the pipeline completes.
- `manifest.path` (string): Relative path under `io.out_root` (default `logs/run_manifest.json`).
