test_that("spin_correlation works for two-map input", {
  set.seed(1)
  x <- rnorm(12)
  y <- x + rnorm(12, sd = 0.2)

  perm_id <- replicate(40, sample.int(12), simplify = "matrix")
  res <- spin_correlation(x, y, perm_id = perm_id, method = "pearson")

  expect_s3_class(res, "spin_correlation_result")
  expect_equal(dim(res$observed), c(1, 1))
  expect_equal(dim(res$p_value), c(1, 1))
  expect_true(is.finite(res$observed[1, 1]))
  expect_true(res$p_value[1, 1] >= 0 && res$p_value[1, 1] <= 1)
})

test_that("spin_correlation works for multi-map input", {
  set.seed(2)
  x <- matrix(rnorm(68 * 2), nrow = 68, ncol = 2, dimnames = list(NULL, c("x1", "x2")))
  y <- matrix(rnorm(68 * 3), nrow = 68, ncol = 3, dimnames = list(NULL, c("y1", "y2", "y3")))
  perm_id <- replicate(30, sample.int(68), simplify = "matrix")

  res <- spin_correlation(x, y, perm_id = perm_id, method = "spearman")
  expect_equal(dim(res$observed), c(2, 3))
  expect_equal(dim(res$p_value), c(2, 3))
  expect_equal(dim(res$null_xy), c(30, 2, 3))
  expect_equal(dim(res$null_yx), c(30, 2, 3))
})

test_that("spin_generate_permutations is reproducible with seed", {
  p1 <- spin_generate_permutations(n_perm = 20, seed = 123)
  p2 <- spin_generate_permutations(n_perm = 20, seed = 123)
  expect_identical(p1, p2)
})

test_that("shape mismatch throws informative error", {
  x <- rnorm(20)
  y <- rnorm(22)
  perm_id <- replicate(10, sample.int(20), simplify = "matrix")
  expect_error(
    spin_correlation(x, y, perm_id = perm_id),
    "E_SPIN_SHAPE"
  )
})
