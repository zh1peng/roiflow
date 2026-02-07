# Spin-correlation examples for roiflow

library(roiflow)

set.seed(42)

# 1) Two-map run ---------------------------------------------------------------
x <- rnorm(68)
y <- 0.4 * x + rnorm(68, sd = 0.8)

res_two <- spin_correlation(
  x = x,
  y = y,
  n_perm = 200,
  method = "pearson",
  seed = 2026
)
print(res_two)
print(res_two$summary)

# 2) Multi-map run -------------------------------------------------------------
x_multi <- cbind(map_a = x, map_b = rnorm(68))
y_multi <- cbind(map_c = y, map_d = rnorm(68), map_e = rnorm(68))

res_multi <- spin_correlation(
  x = x_multi,
  y = y_multi,
  n_perm = 200,
  method = "spearman",
  seed = 2026
)
print(head(spin_results_table(res_multi)))

# 3) Optional saving -----------------------------------------------------------
tmp_out <- file.path(tempdir(), "spin_demo")
save_spin_results(res_multi, out_dir = tmp_out, prefix = "demo", formats = c("csv", "rds"))

# 4) Optional plotting ---------------------------------------------------------
p_heat <- plot_spin_matrix(res_multi, value = "observed")
p_null <- plot_spin_null(res_two)
print(p_heat)
print(p_null)
