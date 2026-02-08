---
name: roiflow-orchestrator
description: Execute `roiflow::run_pipeline()` jobs with strict config contracts. Use when a request involves building, validating, or running a roiflow JSON job file, explaining job JSON parameters, troubleshooting pipeline step failures, or producing deterministic outputs from config-first workflows.
---

# roiflow Orchestrator

Use this skill for config-first execution only (`run_pipeline`/`roiflow`).

## Do Not Hallucinate

Use only this function allowlist unless the user explicitly requests otherwise:

- `read_cfg`
- `validate_cfg`
- `resolve_cfg`
- `run_pipeline`
- `roiflow`

Do not call PET APIs in this skill.

## Preflight Checklist

Run before execution:

1. Verify R executable is available:
   - `Rscript --version` must succeed.
   - If unavailable, stop with: `Rscript not found on PATH`.
2. Verify package availability in current R environment:
   - `requireNamespace("roiflow", quietly = TRUE)` must be `TRUE`, or load from source with `devtools::load_all(".")`.
   - If unavailable, stop with: `Package 'roiflow' is not installed/loaded in this R environment.`
3. Check config path exists (if using file input).
4. Validate config shape (`validate_cfg`) before running.
5. Ensure input CSV path exists (`cfg$io$input`).
6. Ensure ROI detection is feasible:
   - if `roi.cols` is set, all listed columns must exist;
   - if `roi.cols` is null, `roi.regex` must match at least one column.
7. Fail fast with explicit error; do not invent missing columns/steps.

## Environment Bootstrap (if missing)

If preflight fails for runtime setup, provide explicit remediation:

1. Install R and ensure `Rscript` is on PATH.
2. Install dependencies (for example `devtools`).
3. Install or load `roiflow`:
   - GitHub install (if available):
     - `devtools::install_github("<github-owner>/roiflow")`
   - Local source (inside repo):
     - `devtools::load_all(".")`

Recommended runtime guard snippet:

```r
if (Sys.which("Rscript") == "") stop("Rscript not found on PATH", call. = FALSE)
if (!requireNamespace("roiflow", quietly = TRUE)) {
  stop("Package 'roiflow' is not installed in this R library.", call. = FALSE)
}
```

## Fixed Workflow Order (Pipeline Task)

1. Identify job source:
   - existing JSON path, or
   - build new JSON from template.
2. Validate job:
   - `read_cfg` -> `validate_cfg`.
3. Execute:
   - `run_pipeline(cfg_path_or_list)`.
4. Summarize outputs:
   - `out_root`,
   - executed step sequence,
   - key result table row counts,
   - manifest/log paths.

## Job JSON Guidance

Use `references/run-pipeline-job-json.md` for full schema and step parameters.

Use `references/roiflow-contracts.md` for strict input/output contracts and expected output file names.

## Failure Rules

- Invalid step names: stop with explicit unsupported-step error.
- Missing required fields (`io.input`, `pipeline`, etc.): stop with explicit config error.
- Missing input file or required columns: stop with exact missing path/column names.
- Do not auto-correct user intent silently (for example, changing `compare_groups` to `pairwise_groups` without user approval).

## Canonical Examples

### 1) Run an existing job JSON

```r
library(roiflow)
ctx <- run_pipeline("inst/examples/example_cfg.json")
ctx$cfg$io$out_root
ctx$timings
```

### 2) Validate then run

```r
library(roiflow)
cfg <- read_cfg("inst/examples/example_assoc_cfg.json")
cfg <- validate_cfg(cfg)
ctx <- run_pipeline(cfg)
```

### 3) Build a new job JSON from template and run

```r
library(jsonlite)
job <- read_json("inst/examples/job_template.json", simplifyVector = TRUE)
job$io$out_root <- "outputs/my_job"
write_json(job, "tmp_my_job.json", auto_unbox = TRUE, pretty = TRUE)

library(roiflow)
ctx <- run_pipeline("tmp_my_job.json")
```
