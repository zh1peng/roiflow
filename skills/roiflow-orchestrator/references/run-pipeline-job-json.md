# run_pipeline Job JSON Guide

Use this guide when building or reviewing a `run_pipeline` job JSON.

Reference template:

- `inst/examples/job_template.json`
- `inst/examples/job_template.md`

## 1) Top-level Structure

```json
{
  "project": { "name": "my_project", "version": "0.1.0" },
  "io": { ... },
  "columns": { ... },
  "roi": { ... },
  "pipeline": [ ... ],
  "manifest": { ... }
}
```

## 2) `io` Section

- `input` (string): input CSV path.
- `out_root` (string): output root directory.
- `overwrite` (bool): whether existing outputs can be replaced.
- `save_data` (bool): write final dataset.
- `save_tables` (bool): write tables.
- `save_plots` (bool): write plots.
- `plot_formats` (array of strings): e.g., `["png"]`.

## 3) `columns` Section

- `site` (string): ComBat batch column.
- `group` (string): group column or `"auto"`.
- `age`, `sex`, `icv` (string): common covariate mappings.

## 4) `roi` Section

- `regex` (string): ROI detection regex.
- `drop_suffix` (string): optional label cleanup suffix.
- `cols` (array or null): explicit ROI column list; if null, use regex.

## 5) `pipeline` Section

Ordered array of steps:

```json
{ "step": "prep", "params": { ... } }
```

Supported step names:

- `prep`
- `combat`
- `compare_groups`
- `pairwise_groups`
- `assoc_scan`
- `corr_scan`
- `plots`
- `export`

### `prep` params

- `group_ref` (string)
- `sex_ref` (string)
- `factor_levels` (object or null)
- `na_action`: `drop_required` | `drop_all` | `keep`
- `outlier_action`: `none` | `flag` | `clip`
- `verbose` (bool)

### `combat` params

- `covariates` (array)
- `missing_fit`: `impute` | `error`
- `single_batch`: `error` | `noop`
- `verbose` (bool)

### `compare_groups` params

- `covariates` (array)
- `ref` (string)
- `target` (string)
- `model_engine`
- `diagnostics`
- `effect_size`
- `effect_size_ci`
- `effect_size_from`
- `p_adjust_method`
- `family_id`
- `family_desc`

### `pairwise_groups` params

- Similar to `compare_groups`, but no single fixed target pair.

### `assoc_scan` params

- `x` (string, required)
- `covariates` (array)
- `model_engine`
- `diagnostics`
- `effect_size`
- `p_adjust_method`
- `family_id`
- `family_desc`
- `by_group` (bool)
- `interaction` (bool)

### `corr_scan` params

- `x` (string, required)
- `method`: `pearson` | `spearman`
- `covariates` (array)

### `plots` params

- `enabled` (array): any of
  - `brain_map`
  - `pvalue_hist`
  - `top_rois`
  - `raincloud_roi`
  - `scatter_assoc`
- plus plot-specific appearance/settings keys.

### `export` params

- `base_name` (string)
- `formats` (array: `csv` and/or `tsv`)
- `style`: `publication` | `minimal`

## 6) `manifest` Section

- `enabled` (bool)
- `path` (string), usually `logs/run_manifest.json`

## 7) Minimal Working Example

```json
{
  "io": {
    "input": "data/sample.csv",
    "out_root": "outputs/my_run",
    "overwrite": true
  },
  "columns": {
    "site": "Site",
    "group": "auto",
    "age": "Age",
    "sex": "Sex",
    "icv": "ICV"
  },
  "roi": {
    "regex": "^(L_|R_)",
    "drop_suffix": "_thickavg"
  },
  "pipeline": [
    { "step": "prep", "params": { "na_action": "drop_required", "outlier_action": "none" } },
    { "step": "combat", "params": { "missing_fit": "impute", "single_batch": "error" } },
    { "step": "compare_groups", "params": { "ref": "Control", "target": "Case", "p_adjust_method": "fdr" } },
    { "step": "plots", "params": { "enabled": ["brain_map", "pvalue_hist"] } },
    { "step": "export", "params": { "base_name": "ct_group_diff", "formats": ["csv"], "style": "publication" } }
  ],
  "manifest": { "enabled": true, "path": "logs/run_manifest.json" }
}
```

## 8) Build Pattern

1. Copy `inst/examples/job_template.json`.
2. Edit `io.input` and `io.out_root`.
3. Keep only needed steps in `pipeline` in desired order.
4. Fill required params for selected steps (for example `compare_groups.ref/target`).
5. Validate and run with:
   - `cfg <- validate_cfg(read_cfg(path))`
   - `ctx <- run_pipeline(path)`
