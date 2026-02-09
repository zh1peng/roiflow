.make_combat_plot_data <- function() {
  set.seed(11)
  n <- 120
  site <- rep(c("S1", "S2", "S3"), each = 40)
  study <- rep(c("StudyA", "StudyB", "StudyC"), each = 40)
  group <- rep(rep(c("Control", "Case"), each = 20), 3)

  site_effect <- ifelse(site == "S1", -0.6, ifelse(site == "S2", 0.5, 0.1))
  group_effect <- ifelse(group == "Case", 0.35, 0)
  banksts_before <- stats::rnorm(n, sd = 0.5) + site_effect + group_effect
  insula_before <- stats::rnorm(n, sd = 0.45) + 0.5 * site_effect + group_effect

  banksts_after <- banksts_before - 0.85 * site_effect + stats::rnorm(n, sd = 0.08)
  insula_after <- insula_before - 0.8 * site_effect + stats::rnorm(n, sd = 0.08)

  before <- data.frame(
    Site = site,
    Study = study,
    Group = group,
    banksts = banksts_before,
    insula = insula_before,
    stringsAsFactors = FALSE
  )
  after <- data.frame(
    Site = site,
    Study = study,
    Group = group,
    banksts = banksts_after,
    insula = insula_after,
    stringsAsFactors = FALSE
  )

  list(before = before, after = after)
}

.to_long <- function(df) {
  rows <- vector("list", 2)
  rows[[1]] <- data.frame(
    Site = df$Site,
    Study = df$Study,
    Group = df$Group,
    roi = "banksts",
    value = df$banksts,
    stringsAsFactors = FALSE
  )
  rows[[2]] <- data.frame(
    Site = df$Site,
    Study = df$Study,
    Group = df$Group,
    roi = "insula",
    value = df$insula,
    stringsAsFactors = FALSE
  )
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

test_that("combat_plot_roi returns ggplot for wide input with default split", {
  skip_if_not_installed("ggplot2")
  dat <- .make_combat_plot_data()

  p <- combat_plot_roi(
    df_before = dat$before,
    df_after = dat$after,
    roi = "banksts",
    site_col = "Site",
    group_col = "Group",
    study_col = "Study",
    split = FALSE
  )

  expect_s3_class(p, "ggplot")
  expect_true(all(c("before", "after") %in% unique(as.character(p$data$stage))))
  expect_true(all(c("Case", "Control") %in% unique(as.character(p$data$group))))
})

test_that("combat_plot_roi supports split output with two ggplot panels", {
  skip_if_not_installed("ggplot2")
  dat <- .make_combat_plot_data()

  p <- combat_plot_roi(
    df_before = dat$before,
    df_after = dat$after,
    roi = "banksts",
    site_col = "Site",
    group_col = "Group",
    split = TRUE
  )

  expect_s3_class(p, "combat_roi_plot_split")
  expect_s3_class(p$before, "ggplot")
  expect_s3_class(p$after, "ggplot")
})

test_that("combat_plot_roi works for long format and supports site/study filters", {
  skip_if_not_installed("ggplot2")
  dat <- .make_combat_plot_data()
  long_before <- .to_long(dat$before)
  long_after <- .to_long(dat$after)

  p <- combat_plot_roi(
    df_before = long_before,
    df_after = long_after,
    roi = "banksts",
    site_col = "Site",
    group_col = "Group",
    roi_col = "roi",
    value_col = "value",
    study_col = "Study",
    site = "S1",
    study = "StudyA"
  )

  expect_s3_class(p, "ggplot")
  expect_true(all(as.character(p$data$site) == "S1"))
  expect_true(all(as.character(p$data$study) == "StudyA"))
})

test_that("combat_plot_roi errors on multiple ROI names", {
  skip_if_not_installed("ggplot2")
  dat <- .make_combat_plot_data()

  expect_error(
    combat_plot_roi(dat$before, dat$after, roi = c("banksts", "insula")),
    "single non-empty ROI name"
  )
})

test_that("combat_plot_roi errors on unknown ROI", {
  skip_if_not_installed("ggplot2")
  dat <- .make_combat_plot_data()

  expect_error(
    combat_plot_roi(dat$before, dat$after, roi = "not_a_real_roi"),
    "E_COMBAT_PLOT_ROI"
  )
})

test_that("combat_plot_roi warns when group levels are not binary", {
  skip_if_not_installed("ggplot2")
  dat <- .make_combat_plot_data()
  dat3 <- dat$before
  dat3$Group[1:6] <- "Other"
  dat3$Group <- factor(dat3$Group)

  expect_warning(
    combat_plot_roi(dat3, dat$after, roi = "banksts"),
    "W_COMBAT_GROUP_LEVELS"
  )
})
