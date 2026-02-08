# roiflow PET Contracts

This file defines strict contracts for PET workflows.

## 1) Compute API Contract

### Allowed compute APIs

- `pet_available()`
- `pet_corr(...)`
- `pet_corr_save(...)`

### `pet_corr` input contract

- `brain`: numeric vector OR matrix/data.frame (`n_roi x n_maps`)
- PET source:
  - explicit `pet` (vector/matrix/data.frame), or
  - `pet_select` (character names from packaged PET maps)
- control args:
  - `n_perm`
  - `method`
  - `seed`
  - optional `perm_id`
  - optional `hemi`, `na_action`, `return_null`

### Alignment contract

- If row names exist:
  - required ROI names must be present;
  - rows are aligned/reordered by expected ROI order.
- If row names absent:
  - row count must match selected hemisphere/centroid expectations.

### Output contract

- Returns class `pet_corr_result`
- Fields:
  - `results`
  - `nulls` (optional)
  - `meta`
  - `inputs`
  - `call`

Minimum `results` columns:

- `brain_map`
- `pet_map`
- `r_obs`
- `p_spin`
- `n_perm`
- `method`
- `seed`

## 2) Plot API Contract

Allowed plotting APIs:

- `pet_plot_scatter(pet_corr_obj, brain_map, pet_map, ...)`
- `pet_plot_bar(pet_corr_obj, ...)`
- `pet_plot_radar(pet_corr_obj, ...)`

Rules:

- Input must be `pet_corr_result`.
- Requested map names must exist in result object.
- Return `ggplot` object.
- No side-effect saving unless explicitly requested.

## 3) Save Contract

Using `pet_corr_save` (or `pet_corr(save=TRUE)`):

- Outputs:
  - `<out_dir>/<out_prefix>.csv`
  - `<out_dir>/<out_prefix>.rds`
- Defaults:
  - `out_prefix = "pet_corr_results"`

## 4) Failure Contract

Do not guess.

- Unknown `pet_select` names -> explicit selection error.
- Shape/alignment mismatch -> explicit mismatch error.
- Missing required objects/data -> explicit missing-object/path error.
- Missing map names for plotting -> explicit map-not-found error.
