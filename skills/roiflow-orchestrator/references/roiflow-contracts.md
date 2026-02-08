# roiflow run_pipeline Contracts

This contract applies to the `roiflow-orchestrator` skill.

## 1) Accepted Execution APIs

- `run_pipeline(cfg_or_path)`
- `roiflow(cfg_or_path)` (alias)
- Optional validation helpers:
  - `read_cfg(path)`
  - `validate_cfg(cfg)`
  - `resolve_cfg(cfg, df)` for debug-only inspection

No PET APIs are in scope for this skill.

## 2) Input Contract

`cfg_or_path` must be either:

1. JSON config path string, or
2. Config list object compatible with `validate_cfg()`.

Required minimum config components:

- `io.input`
- `io.out_root`
- `columns.site` and other column mappings used by selected steps
- `roi.regex` or `roi.cols`
- non-empty `pipeline` step array

Runtime preconditions:

- `Rscript --version` succeeds.
- `requireNamespace("roiflow", quietly = TRUE)` is `TRUE` (or package loaded from source via `devtools::load_all(".")`).

## 3) Output Contract

`run_pipeline` returns context object `ctx` with fields:

- `cfg`
- `input`
- `data`
- `specs`
- `results`
- `plots`
- `exports`
- `manifest`
- `log`
- `timings`
- `paths`
- `meta`

## 4) Saved Output Contract (when pipeline includes `export`)

Files are controlled by config and active steps.

Stable defaults:

- manifest: `<out_root>/logs/run_manifest.json`
- run log: `<out_root>/logs/run_log.json`
- data: `<out_root>/data/...`
- tables: `<out_root>/tables/...`
- figures: `<out_root>/figures/...`

Common table filenames produced by step combinations:

- compare-groups export:
  - `<out_root>/tables/<base_name>_raw.csv|tsv`
  - `<out_root>/tables/<base_name>_formatted.csv|tsv`
  - `<out_root>/tables/ct_group_summary.csv` (if summaries available)
  - `<out_root>/tables/ct_diagnostics.csv` (if diagnostics available)
- pairwise-groups export:
  - `<out_root>/tables/pairwise_groups_raw.csv|tsv`
  - `<out_root>/tables/pairwise_groups_formatted.csv|tsv`
- assoc-scan export:
  - `<out_root>/tables/assoc_scan.csv`
- corr-scan export:
  - `<out_root>/tables/corr_scan.csv`

## 5) Failure Contract

Fail fast; do not guess.

- Invalid step names -> explicit unsupported-step error.
- Missing input file -> explicit file-not-found error.
- Missing columns/ROI matches -> explicit missing-column/ROI error.
- Existing files with `overwrite = false` -> explicit exists-conflict error.
