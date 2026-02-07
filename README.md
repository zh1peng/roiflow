# roiflow

**Workflow for ROI-Based Neuroimaging Data Analysis**

An R package providing a comprehensive, configurable workflow for preprocessing, analyzing, and clustering ROI-based neuroimaging data.

## Installation

```r
# Install from local source
devtools::install_local("path/to/roiflow")

# Or install development version from GitHub (once published)
# devtools::install_github("username/roiflow")
```

## Features

- **Flexible Configuration**: Customize preprocessing via `prep_spec()`
- **Comprehensive Logging**: Track all operations with structured log codes
- **Quality Control Metrics**: Detailed QC statistics for data validation
- **Multiple Outlier Methods**: IQR or SD-based detection with configurable thresholds
- **Smart Missing Data Handling**: Drop required columns only, all columns, or keep all data
- **Site Harmonization**: Automatic reference level setting for multi-site studies
- **LLM-Agent Friendly**: Clear function names, comprehensive docs, informative messages

## Quick Start

### Basic Usage

```r
library(roiflow)

# Load your ROI data
# Expected columns: L_*, R_* (ROI measurements), Age, Sex, Site, Group, ICV
raw_data <- read.csv("your_roi_data.csv")

# Preprocess with default settings - returns cleaned data
clean_data <- prep(raw_data)
# prep(): 1000 -> 987 rows; NA-dropped=13; outliers=142 points
```

### Get Full Results with QC Metrics

```r
# Return full result object with QC metrics and logs
result <- prep(raw_data, return = "result")
print(result)
# <roiflow_prep_result>
# Rows: 1000 -> 987 (dropped NA: 13)
# ROI cols: 68 | Outlier points: 142 | New NA from numeric: 0

# Access components
clean_data <- result$data
qc_metrics <- result$qc
logs <- result$log
spec_used <- result$spec

# Inspect QC metrics
str(result$qc)
# List of 11
#  $ n_rows_in           : int 1000
#  $ n_rows_out          : int 987
#  $ n_roi_cols          : int 68
#  $ n_rows_dropped_na   : int 13
#  $ n_zero_to_na        : int 45
#  $ n_outlier_points    : int 142
#  $ outliers_per_roi    : Named int [1:68] ...
#  $ required_cols       : chr [1:72] ...

# Inspect logs
lapply(result$log, function(x) x$message)
```

### Custom Configuration

```r
# Create custom preprocessing specification
spec <- prep_spec(
  roi_regex = "^(L_|R_)",           # Pattern to identify ROI columns
  outlier_k = 2.5,                   # More conservative outlier detection
  outlier_method = "iqr",            # Use IQR method (default)
  na_action = "drop_required",       # Drop rows with NA in required cols only
  site_ref = "SITE_A",               # Use specific site as reference
  icv_scale = 1e6,                   # Rescale ICV to millions
  factor_levels = list(              # Custom factor reference levels
    Sex = "0",
    Group = "0",
    Diagnosis = "HC"
  )
)

# Use custom spec
clean_data <- prep(raw_data, spec = spec)
```

### Common Use Cases

```r
# 1. Keep all data, just flag outliers (no clipping)
spec <- prep_spec(
  na_action = "keep",
  outlier_action = "flag"
)
result <- prep(raw_data, spec = spec, return = "result")
# Check which ROIs have outliers
result$qc$outliers_per_roi

# 2. More aggressive outlier removal (k=2 instead of 3)
spec <- prep_spec(outlier_k = 2)
clean_data <- prep(raw_data, spec = spec)

# 3. Use SD method instead of IQR for outliers
spec <- prep_spec(outlier_method = "sd", outlier_k = 3)
clean_data <- prep(raw_data, spec = spec)

# 4. Drop rows with ANY missing data
spec <- prep_spec(na_action = "drop_all")
clean_data <- prep(raw_data, spec = spec)

# 5. Strict mode - stop on warnings
clean_data <- prep(raw_data, strict = TRUE)

# 6. Silent mode - no progress messages
clean_data <- prep(raw_data, verbose = FALSE)
```

## Config-First Pipeline (Agent-Ready)

For end-to-end runs (prep -> ComBat -> analysis -> plots -> export), `roiflow` supports
a config-first pipeline driven by a single JSON file.

- Example config: `inst/examples/example_cfg.json`
- Job template + parameter documentation: `inst/examples/job_template.json`, `inst/examples/job_template.md`
- API reference (Usage + Arguments): `inst/examples/api_reference.md`
- Single entrypoint: `run_pipeline("path/to/cfg.json")` (alias: `roiflow("path/to/cfg.json")`)

### Run The Demos (From Project Root)

```bash
Rscript scripts/demo_two_group_human.R
```

```bash
Rscript scripts/demo_two_group_agent.R
```

```bash
Rscript scripts/demo_mutiple_group_human.R
```

```bash
Rscript scripts/demo_mutiple_group_agent.R
```

```bash
Rscript scripts/demo_assoc_human.R
```

```bash
Rscript scripts/demo_assoc_agent.R
```

### JSON Pipeline (One Call)

```r
devtools::load_all(".")
ctx <- run_pipeline("inst/examples/example_cfg.json")
```

Association JSON example:

```r
devtools::load_all(".")
ctx <- run_pipeline("inst/examples/example_assoc_cfg.json")
```

Outputs are written under `outputs/<project>/` with subfolders:

- `data/` harmonized datasets
- `tables/` raw + formatted results tables
- `figures/` plots (PNG by default)
- `logs/` run manifest + structured log

### Results Table Formatting

The pipeline uses package-level utilities (not ad-hoc demo code):

- `format_results_table(res_tbl, style = c("minimal","publication"), ...)`
- `export_results_tables(res_tbl, out_dir, base_name = "...", ...)`

## Plotting Utilities

All plot utilities live in `R/plot.R` and return `ggplot` objects (saving is
handled by demos or the pipeline export helpers).

| Function | What It Plots | Required Inputs (Typical) |
|---|---|---|
| `plot_dot()` | group means ± SE (dot + error bars) | `summary_tbl` with `var`, `groups`, `mean`, `se` |
| `plot_bar()` | group means ± SE (bar + error bars) | `summary_tbl` with `var`, `groups`, `mean`, `se` |
| `plot_violin()` | raw distributions (violin + box) | `data`, `x`, `y` |
| `plot_radar()` | per-group radar profiles | `summary_tbl` with `var`, `groups`, `mean` |
| `plot_quadrant()` | x vs y scatter with quadrants | `stat_tbl`, `x_col`, `y_col` |
| `plot_brain_map()` | ggseg brain map from atlas labels | `tbl`, `atlas_spec` (`atlas`, `label_col`, `value_col`) |
| `plot_brain_map_results()` | brain map from ROI names (L_/R_ -> ggseg labels) | `res_tbl` with `var` and a `value_col` |
| `plot_qc_pvalue_hist()` | p-value histogram | `res_tbl` with `p_col` |
| `plot_top_rois()` | top-N ROIs by |metric| | `res_tbl` with `metric` and label column |
| `plot_auc_hist()` | AUC null/true histograms | numeric vectors `auc_null`, optional `auc_true` |
| `plot_raincloud_roi()` | ROI distribution (2+ groups; raw or marginal; optional pairwise markers) | `data`, `roi`, `group_col` |
| `plot_scatter_assoc()` | ROI vs predictor scatter (raw or marginal) | `data`, `roi`, `x` |

## PET Spin-Correlation Utilities

`roiflow` also includes PET DK-68 maps and spin-test helpers:

- `pet_available()` to list packaged PET targets + metadata
- `pet_corr()` to correlate 1..N brain maps vs 1..N PET maps with `p_spin`
- `pet_plot_scatter()`, `pet_plot_bar()`, `pet_plot_radar()` for publication-ready figures

Quick example:

```r
data("pet_maps_dk68")

# Example brain maps (replace with your own DK-68 ROI vectors/matrix)
brain <- pet_maps_dk68[, c("D1", "DAT")]

# Fast demo permutation matrix (use larger n_perm for analysis)
perm_id <- replicate(100, sample.int(nrow(brain)), simplify = "matrix")

res <- pet_corr(
  brain = brain,
  pet_select = c("D1", "D2", "DAT", "MOR"),
  perm_id = perm_id,
  n_perm = 100
)

head(res$results)
pet_plot_scatter(res, brain_map = "D1", pet_map = "D1")
pet_plot_bar(res, style = "bar")
pet_plot_radar(res, metric = "r_obs")
```

## Main Functions

### `prep()`

Main preprocessing function that performs a complete pipeline:

1. **Column Identification**: Detects ROI, numeric, and zero-to-NA columns via regex
2. **Type Conversion**: Converts columns to numeric with NA tracking
3. **ICV Rescaling**: Divides ICV by scale factor (default: 1,000,000)
4. **Factor Handling**: Sets categorical variables with appropriate reference levels
5. **Zero to NA**: Converts zeros to NA in ROI/ICV columns
6. **Missing Data**: Handles NA according to policy (drop required/all/keep)
7. **Outlier Detection**: Clips or flags outliers using IQR or SD method

**Parameters:**
- `df`: Data frame with ROI measurements
- `spec`: Preprocessing specification from `prep_spec()`
- `return`: `"data"` (default) or `"result"` (includes QC and logs)
- `verbose`: Print progress messages (default: `TRUE` in interactive sessions)
- `strict`: Stop on warnings (default: `FALSE`)

### `prep_spec()`

Creates a configuration object controlling all preprocessing behavior.

**Key Parameters:**
- `roi_regex`: Pattern to identify ROI columns (default: `"^(L_|R_)"`)
- `outlier_k`: Multiplier for outlier fences (default: `3`)
- `outlier_method`: `"iqr"` or `"sd"` (default: `"iqr"`)
- `outlier_action`: `"clip"`, `"flag"`, or `"none"` (default: `"clip"`)
- `na_action`: `"drop_required"`, `"drop_all"`, or `"keep"` (default: `"drop_required"`)
- `site_ref`: `"largest"`, `"keep"`, or specific site name (default: `"largest"`)
- `icv_scale`: Factor to divide ICV by (default: `1e6`)
- `factor_levels`: Named list of factor reference levels

## Expected Data Format

Your input data frame should contain:

- **ROI columns**: Named with `L_*` or `R_*` prefix (e.g., `L_hippocampus`, `R_amygdala`)
- **Age**: Participant age (optional)
- **Sex**: Binary sex variable, coded as 0/1 (optional)
- **Site**: Scanning site identifier (optional)
- **Group**: Group membership, coded as 0/1 (optional)
- **ICV**: Intracranial volume (optional)

## Quality Control Metrics

When using `return = "result"`, you get comprehensive QC metrics:

```r
result <- prep(raw_data, return = "result")

# QC metrics available:
result$qc$n_rows_in              # Input row count
result$qc$n_rows_out             # Output row count
result$qc$n_cols_in              # Input column count
result$qc$n_cols_out             # Output column count
result$qc$n_roi_cols             # Number of ROI columns detected
result$qc$n_rows_dropped_na      # Rows removed due to missing data
result$qc$n_zero_to_na           # Zeros converted to NA
result$qc$n_outlier_points       # Total outlier values detected/clipped
result$qc$outliers_per_roi       # Named vector: outliers per ROI
result$qc$required_cols          # Columns used for completeness check
result$qc$n_na_new_from_numeric  # NAs introduced by type coercion
```

## Structured Logging

The package uses structured log codes for programmatic parsing:

```r
result <- prep(raw_data, return = "result")

# Log structure
result$log[[1]]
# $code
# [1] "I_ICV_RESCALE"
#
# $message
# [1] "Rescaled ICV by /1e+06."
#
# $context
# $context$icv_col
# [1] "ICV"
# $context$icv_scale
# [1] 1e+06

# Log codes:
# W_* = Warnings (e.g., W_NO_ROI_COLS, W_NUMERIC_COERCE_NA, W_FACTOR_REF_MISSING)
# I_* = Informational (e.g., I_ICV_RESCALE, I_FACTOR_RELEVEL, I_SITE_REF_LARGEST)
```

## LLM Agent Friendly Design

This package is designed to be easily understood and used by both humans and LLM agents:

- **Clear function names**: `prep()`, `prep_spec()` - concise and descriptive
- **Comprehensive documentation**: Every function has detailed roxygen2 docs
- **Flexible configuration**: Separate spec object for easy customization
- **Structured logging**: Machine-parseable log codes with context
- **QC metrics**: Quantitative validation of preprocessing steps
- **Informative messages**: Verbose mode provides step-by-step progress
- **Input validation**: Clear error messages for invalid inputs
- **Sensible defaults**: Works out-of-the-box for common use cases
- **Return options**: Get just data or full result object with metadata

## Development

### Building Documentation

```r
# Generate documentation from roxygen2 comments
devtools::document()
```

### Running Tests

```r
# Run package tests (once tests are added)
devtools::test()
```

### Checking Package

```r
# Run R CMD check
devtools::check()
```

### Loading Package During Development

```r
# Load package from source
devtools::load_all()

# Or install locally
devtools::install()
```

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.

## License

MIT License - see LICENSE file for details.

## Citation

If you use this package in your research, please cite:

```
[Citation information to be added]
```

## Contact

[Contact information to be added]
