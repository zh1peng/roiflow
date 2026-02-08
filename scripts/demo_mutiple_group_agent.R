#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("devtools is required to run this demo from source. Install with install.packages('devtools').", call. = FALSE)
  }
  devtools::load_all(".")
})

# Agent demo (multi-group/pairwise, config-first):
# Reads one job JSON and runs the full pipeline via run_pipeline().
#
# Usage:
#   Rscript scripts/demo_mutiple_group_agent.R
#   Rscript scripts/demo_mutiple_group_agent.R inst/examples/example_multigroup_cfg.json

args <- commandArgs(trailingOnly = TRUE)
cfg_path <- if (length(args) >= 1 && nzchar(args[1])) {
  args[1]
} else {
  file.path("inst", "examples", "example_multigroup_cfg.json")
}

if (!file.exists(cfg_path)) {
  stop(sprintf("Config file not found: %s", cfg_path), call. = FALSE)
}

t0 <- Sys.time()
cat("\n=== Multi-Group Demo (Agent, Config-First) ===\n")
cat(sprintf("Config: %s\n", cfg_path))

ctx <- run_pipeline(cfg_path)

t1 <- Sys.time()
elapsed <- as.numeric(difftime(t1, t0, units = "secs"))

cat(sprintf("Output root: %s\n", ctx$cfg$io$out_root))
cat(sprintf("Pipeline steps: %s\n", paste(vapply(ctx$cfg$pipeline, `[[`, character(1), "step"), collapse = " -> ")))
rows_in <- if (!is.null(ctx$meta$n_rows_input)) ctx$meta$n_rows_input else NA_integer_
cat(sprintf("Rows: %d -> %d\n", rows_in, nrow(ctx$data)))

group_col <- if (!is.null(ctx$cfg$resolved$group_col)) {
  ctx$cfg$resolved$group_col
} else {
  ctx$cfg$columns$group
}
if (!is.null(group_col) && group_col %in% names(ctx$data)) {
  cat("Group counts:\n")
  print(table(ctx$data[[group_col]]))
}

if (!is.null(ctx$results$pairwise_groups$data)) {
  res_tbl <- ctx$results$pairwise_groups$data
  cat(sprintf("Pairwise rows: %d\n", nrow(res_tbl)))
  if ("p_adj" %in% names(res_tbl)) {
    cat(sprintf("Min FDR p: %.4g\n", min(res_tbl$p_adj, na.rm = TRUE)))
  }
} else if (!is.null(ctx$results$compare_groups$data)) {
  res_tbl <- ctx$results$compare_groups$data
  cat(sprintf("Compare groups rows: %d\n", nrow(res_tbl)))
  if ("p_adj" %in% names(res_tbl)) {
    cat(sprintf("Min FDR p: %.4g\n", min(res_tbl$p_adj, na.rm = TRUE)))
  }
}

manifest_rel <- if (!is.null(ctx$cfg$manifest$path)) ctx$cfg$manifest$path else file.path("logs", "run_manifest.json")
manifest_path <- file.path(ctx$cfg$io$out_root, manifest_rel)
log_path <- file.path(ctx$cfg$io$out_root, "logs", "run_log.json")
cat(sprintf("Manifest: %s\n", manifest_path))
cat(sprintf("Run log: %s\n", log_path))
cat(sprintf("Total seconds: %.2f\n", elapsed))
