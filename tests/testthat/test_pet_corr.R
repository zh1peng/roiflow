test_that("pet_corr single brain map vs single PET map returns expected columns", {
  data("pet_maps_dk68")

  brain <- pet_maps_dk68[, "D1"] + stats::rnorm(nrow(pet_maps_dk68), sd = 0.01)
  perm_id <- replicate(20, sample.int(nrow(pet_maps_dk68)), simplify = "matrix")

  res <- pet_corr(
    brain = brain,
    pet_select = "D1",
    perm_id = perm_id,
    n_perm = 20,
    return_null = TRUE
  )

  expect_s3_class(res, "pet_corr_result")
  expect_equal(nrow(res$results), 1L)
  expect_true(all(c(
    "brain_map", "pet_map", "r_obs", "p_spin", "n_perm", "method", "seed"
  ) %in% names(res$results)))
  expect_true(is.data.frame(res$nulls))
  expect_equal(nrow(res$nulls), 1L)
})

test_that("pet_corr multi brain map vs multi PET maps returns cartesian rows", {
  data("pet_maps_dk68")

  brain <- pet_maps_dk68[, c("D1", "MOR"), drop = FALSE]
  pet <- pet_maps_dk68[, c("D1", "DAT", "MOR"), drop = FALSE]
  perm_id <- replicate(15, sample.int(nrow(pet_maps_dk68)), simplify = "matrix")

  res <- pet_corr(
    brain = brain,
    pet = pet,
    perm_id = perm_id,
    n_perm = 15,
    return_null = FALSE
  )

  expect_equal(nrow(res$results), ncol(brain) * ncol(pet))
  expect_setequal(unique(res$results$brain_map), colnames(brain))
  expect_setequal(unique(res$results$pet_map), colnames(pet))
})

test_that("pet_corr spin p-values are reproducible with the same seed", {
  data("pet_maps_dk68")

  brain <- pet_maps_dk68[, c("D1", "DAT"), drop = FALSE]
  set.seed(99)
  res1 <- pet_corr(
    brain = brain,
    pet_select = c("D1", "D2", "DAT"),
    n_perm = 20,
    seed = 1234,
    return_null = FALSE
  )
  set.seed(99)
  res2 <- pet_corr(
    brain = brain,
    pet_select = c("D1", "D2", "DAT"),
    n_perm = 20,
    seed = 1234,
    return_null = FALSE
  )

  expect_equal(res1$results$p_spin, res2$results$p_spin)
})

test_that("pet_corr errors on invalid pet_select names", {
  data("pet_maps_dk68")

  expect_error(
    pet_corr(
      brain = pet_maps_dk68[, "D1"],
      pet_select = "NOT_A_REAL_PET_MAP",
      n_perm = 10
    ),
    "E_PET_SELECT"
  )
})

test_that("pet_available returns packaged PET metadata", {
  tbl <- pet_available()
  expect_true(is.data.frame(tbl))
  expect_true(all(c("pet_map", "target", "modality", "parcellation") %in% names(tbl)))
  expect_true(nrow(tbl) >= 1)
})
