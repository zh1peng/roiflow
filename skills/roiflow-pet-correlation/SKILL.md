---
name: roiflow-pet-correlation
description: Run roiflow PET map workflows with strict validation. Use when requests involve listing packaged PET maps, running `pet_corr` for brain map vs PET map spin-correlation, saving PET correlation outputs, or plotting PET correlation results with `pet_plot_scatter`, `pet_plot_bar`, or `pet_plot_radar`.
---

# roiflow PET Correlation

Use this skill for PET map correlation and PET-result plotting only.

## Do Not Hallucinate

Use only this allowlist unless user explicitly requests otherwise:

- `pet_available`
- `pet_corr`
- `pet_corr_save`
- `pet_plot_scatter`
- `pet_plot_bar`
- `pet_plot_radar`

Do not use `run_pipeline` in this skill.

## Preflight Checklist

1. Load package (`library(roiflow)` or `devtools::load_all(".")`).
2. Validate inputs:
   - `brain` is numeric vector/matrix/data.frame.
   - `pet` or `pet_select` is provided according to task.
3. Validate PET map selection:
   - if using `pet_select`, names must exist in `pet_available()`.
4. Validate alignment:
   - matching ROI rows/order between `brain` and PET maps.
   - if row names exist, enforce expected ROI names; if absent, enforce row-count compatibility.
5. Fail fast on errors; do not guess names/ordering.

Read strict contracts in `references/pet-contracts.md`.

## Fixed Workflow Order

### Task A: PET correlation compute

1. Preflight checklist.
2. Resolve PET maps:
   - explicit `pet`, or
   - `pet_select` against packaged PET metadata.
3. Run `pet_corr(...)` with explicit `n_perm`, `method`, and `seed` when reproducibility matters.
4. Return and summarize `results`.
5. Save only if user requests (`save=TRUE` or `pet_corr_save()`).

### Task B: PET plotting

1. Require valid `pet_corr_result`.
2. Confirm requested map names exist.
3. Plot with exactly one of:
   - `pet_plot_scatter`
   - `pet_plot_bar`
   - `pet_plot_radar`
4. Return plot object; save only if explicitly requested.

## Failure Rules

- Invalid `pet_select`: stop with explicit selection error, do not guess replacements.
- Invalid shape/alignment: stop with explicit mismatch error, no truncation/padding.
- Missing map name in plotting call: stop with explicit map-not-found error.

## Canonical Examples

### 1) PET map listing

```r
library(roiflow)
pet_available()
```

### 2) Correlation run

```r
library(roiflow)
data("pet_maps_dk68")
perm_id <- replicate(100, sample.int(nrow(pet_maps_dk68)), simplify = "matrix")
res <- pet_corr(
  brain = pet_maps_dk68[, c("D1", "DAT")],
  pet_select = c("D1", "D2", "DAT"),
  perm_id = perm_id,
  n_perm = 100,
  method = "pearson",
  seed = 123
)
head(res$results)
```

### 3) Plot and save

```r
p <- pet_plot_scatter(res, brain_map = "D1", pet_map = "D2")
files <- pet_corr_save(res, out_dir = "outputs/pet_demo", out_prefix = "pet_corr_results")
files
```
