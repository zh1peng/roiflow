# Plotting utilities (consume tables, no model fitting)

.as_plot_spec <- function(plot_spec) {
  if (is.null(plot_spec)) return(list())
  if (!is.list(plot_spec)) stop("plot_spec must be a list.", call. = FALSE)
  plot_spec
}

.apply_common_plot_spec <- function(p, plot_spec) {
  plot_spec <- .as_plot_spec(plot_spec)
  if (length(plot_spec) == 0) return(p)

  .check_pkg("ggplot2", "plotting")

  if (!is.null(plot_spec$theme)) {
    base_size <- plot_spec$base_size %||% 11
    th <- switch(
      plot_spec$theme,
      minimal = ggplot2::theme_minimal(base_size = base_size),
      classic = ggplot2::theme_classic(base_size = base_size),
      bw = ggplot2::theme_bw(base_size = base_size),
      none = ggplot2::theme(),
      stop("Unknown plot_spec$theme (use: minimal/classic/bw/none).", call. = FALSE)
    )
    p <- p + th
  } else if (!is.null(plot_spec$base_size)) {
    p <- p + ggplot2::theme(text = ggplot2::element_text(size = plot_spec$base_size))
  }

  if (!is.null(plot_spec$legend_position)) {
    p <- p + ggplot2::theme(legend.position = plot_spec$legend_position)
  }
  if (!is.null(plot_spec$title)) p <- p + ggplot2::labs(title = plot_spec$title)
  if (!is.null(plot_spec$subtitle)) p <- p + ggplot2::labs(subtitle = plot_spec$subtitle)
  if (!is.null(plot_spec$caption)) p <- p + ggplot2::labs(caption = plot_spec$caption)
  if (!is.null(plot_spec$xlab)) p <- p + ggplot2::labs(x = plot_spec$xlab)
  if (!is.null(plot_spec$ylab)) p <- p + ggplot2::labs(y = plot_spec$ylab)

  p
}

.apply_discrete_scales <- function(p, plot_spec, use_color = FALSE, use_fill = FALSE) {
  plot_spec <- .as_plot_spec(plot_spec)
  if (length(plot_spec) == 0) return(p)
  .check_pkg("ggplot2", "plotting")

  if (use_color && !is.null(plot_spec$colors)) {
    p <- p + ggplot2::scale_color_manual(values = plot_spec$colors)
  }
  if (use_fill) {
    vals <- plot_spec$fill_colors %||% plot_spec$colors %||% NULL
    if (!is.null(vals)) {
      p <- p + ggplot2::scale_fill_manual(values = vals)
    }
  }
  p
}

.roi_to_ggseg_label <- function(roi, drop_suffix = "_thickavg") {
  lab <- as.character(roi)
  lab <- sub("^L_", "lh_", lab)
  lab <- sub("^R_", "rh_", lab)
  if (!is.null(drop_suffix) && nzchar(drop_suffix)) {
    lab <- sub(paste0(drop_suffix, "$"), "", lab)
  }
  lab
}

.fit_lm_complete <- function(df, formula, cols_needed) {
  cols_needed <- unique(cols_needed)
  cols_needed <- cols_needed[!is.na(cols_needed) & nzchar(cols_needed)]
  .ensure_cols(df, cols_needed, context = "plot")
  keep <- stats::complete.cases(df[, cols_needed, drop = FALSE])
  d <- df[keep, , drop = FALSE]
  if (nrow(d) == 0) stop("No complete cases for plotting.", call. = FALSE)
  fit <- stats::lm(formula, data = d)
  list(fit = fit, data = d, keep = keep)
}

.adjust_outcome_keep_terms <- function(fit, data, y_col, keep_terms, center = TRUE) {
  y <- data[[y_col]]
  mm <- stats::model.matrix(fit)
  beta <- stats::coef(fit)
  beta[is.na(beta)] <- 0

  trm <- stats::terms(fit)
  term_labels <- attr(trm, "term.labels") %||% character(0)
  assign <- attr(mm, "assign")

  keep_terms <- keep_terms %||% character(0)
  keep_idx <- match(keep_terms, term_labels)
  keep_idx <- keep_idx[!is.na(keep_idx)]
  remove_idx <- setdiff(seq_along(term_labels), keep_idx)
  remove_cols <- which(assign %in% remove_idx)

  eff_remove <- 0
  if (length(remove_cols) > 0) {
    eff_remove <- as.numeric(mm[, remove_cols, drop = FALSE] %*% beta[remove_cols])
  }

  y_adj <- y - eff_remove
  if (isTRUE(center)) {
    y_adj <- y_adj + mean(eff_remove, na.rm = TRUE)
  }
  y_adj
}

# Convert p-values to significance stars (default cutpoints).
.p_to_stars <- function(p) {
  if (is.na(p) || !is.finite(p)) return(NA_character_)
  if (p <= 0.001) return("***")
  if (p <= 0.01) return("**")
  if (p <= 0.05) return("*")
  if (p <= 0.1) return(".")
  "ns"
}

#' Pairwise Tests for a Single ROI Across Groups
#'
#' @description
#' Convenience utility used by the plotting layer to compute pairwise
#' comparisons for a single ROI across all group levels.
#'
#' Supports simple unadjusted tests (\code{t_test}/\code{wilcox}) and an
#' \code{lm}-based test that adjusts for covariates by fitting
#' \code{ROI ~ group + covariates} on each pair (target - ref).
#'
#' @param data Data frame.
#' @param roi Character. ROI column name.
#' @param group_col Character. Grouping column (>= 2 levels).
#' @param covariates Optional character vector of covariate column names.
#' @param value Character. \code{"raw"} or \code{"marginal"}. Only used for
#'   \code{t_test}/\code{wilcox} (marginal values are derived from
#'   \code{ROI ~ group + covariates} by removing covariate contributions while
#'   keeping group effects).
#' @param method Character. One of \code{"t_test"}, \code{"wilcox"}, \code{"lm"}.
#' @param p_adjust Character. Passed to \code{\link[stats]{p.adjust}}. Use
#'   \code{"none"} to skip adjustment.
#' @param pairs Optional list of length-2 character vectors specifying which
#'   pairs to test. Default: all pairs.
#'
#' @return A data frame with columns: \code{group1}, \code{group2},
#'   \code{n1}, \code{n2}, \code{estimate}, \code{statistic},
#'   \code{p_value}, \code{p_adj}, \code{method}.
#' @export
pairwise_tests_roi <- function(data, roi, group_col,
                               covariates = NULL,
                               value = c("raw", "marginal"),
                               method = c("t_test", "wilcox", "lm"),
                               p_adjust = c("holm", "fdr", "bonferroni", "none"),
                               pairs = NULL) {
  value <- match.arg(value)
  method <- match.arg(method)
  p_adjust <- match.arg(p_adjust)

  stopifnot(is.data.frame(data))
  .ensure_cols(data, c(roi, group_col, covariates), context = "pairwise_tests_roi")

  g0 <- data[[group_col]]
  g_levels <- if (is.factor(g0)) levels(droplevels(g0)) else unique(as.character(g0))
  if (length(g_levels) < 2) stop("pairwise_tests_roi requires >= 2 group levels.", call. = FALSE)

  if (is.null(pairs)) {
    pairs <- utils::combn(g_levels, 2, simplify = FALSE)
  } else {
    if (!is.list(pairs) || any(vapply(pairs, length, integer(1)) != 2)) {
      stop("pairs must be a list of length-2 character vectors.", call. = FALSE)
    }
  }

  # Precompute marginal values if requested (for non-model tests only).
  y_use <- NULL
  d_y <- NULL
  if (method %in% c("t_test", "wilcox") && identical(value, "marginal") &&
      !is.null(covariates) && length(covariates) > 0) {
    cols_needed <- unique(c(roi, group_col, covariates))
    keep <- stats::complete.cases(data[, cols_needed, drop = FALSE])
    d_y <- data[keep, , drop = FALSE]
    d_y[[group_col]] <- droplevels(factor(as.character(d_y[[group_col]]), levels = g_levels))
    rhs <- paste(c(.quote_var(group_col), vapply(covariates, .quote_var, character(1))), collapse = " + ")
    f <- stats::as.formula(paste(.quote_var(roi), "~", rhs))
    fit <- stats::lm(f, data = d_y)
    y_use <- .adjust_outcome_keep_terms(fit, d_y, y_col = roi, keep_terms = group_col, center = TRUE)
  }

  out <- lapply(pairs, function(pair) {
    ref <- pair[1]
    target <- pair[2]

    cols_needed <- unique(c(roi, group_col, covariates))
    d <- data[data[[group_col]] %in% c(ref, target), , drop = FALSE]
    d[[group_col]] <- factor(as.character(d[[group_col]]), levels = c(ref, target))

    keep <- stats::complete.cases(d[, cols_needed, drop = FALSE])
    d <- d[keep, , drop = FALSE]
    n1 <- sum(d[[group_col]] == ref)
    n2 <- sum(d[[group_col]] == target)
    if (n1 == 0 || n2 == 0) {
      return(data.frame(
        group1 = ref, group2 = target,
        n1 = n1, n2 = n2,
        estimate = NA_real_, statistic = NA_real_,
        p_value = NA_real_,
        method = method,
        stringsAsFactors = FALSE
      ))
    }

    if (method %in% c("t_test", "wilcox")) {
      if (!is.null(y_use) && !is.null(d_y)) {
        d2 <- d_y[d_y[[group_col]] %in% c(ref, target), , drop = FALSE]
        d2[[group_col]] <- factor(as.character(d2[[group_col]]), levels = c(ref, target))
        y <- y_use[d_y[[group_col]] %in% c(ref, target)]
      } else {
        d2 <- d
        y <- d2[[roi]]
      }
      g <- d2[[group_col]]
      if (method == "t_test") {
        tt <- stats::t.test(y ~ g)
        est <- unname(diff(tt$estimate)) # mean(target) - mean(ref)
        stat <- unname(tt$statistic)
        pv <- tt$p.value
      } else {
        wt <- stats::wilcox.test(y ~ g, exact = FALSE)
        est <- unname(stats::median(y[g == target], na.rm = TRUE) - stats::median(y[g == ref], na.rm = TRUE))
        stat <- unname(wt$statistic)
        pv <- wt$p.value
      }
      return(data.frame(
        group1 = ref, group2 = target,
        n1 = sum(g == ref), n2 = sum(g == target),
        estimate = est,
        statistic = stat,
        p_value = pv,
        method = method,
        stringsAsFactors = FALSE
      ))
    }

    # method == "lm"
    rhs <- paste(c(.quote_var(group_col), vapply(covariates %||% character(0), .quote_var, character(1))), collapse = " + ")
    f <- stats::as.formula(paste(.quote_var(roi), "~", rhs))
    fit <- stats::lm(f, data = d)
    sm <- summary(fit)
    coef_names <- rownames(sm$coefficients)
    idx <- grep(paste0("^", group_col), coef_names)
    idx <- setdiff(idx, which(coef_names == "(Intercept)"))
    if (length(idx) != 1) {
      return(data.frame(
        group1 = ref, group2 = target,
        n1 = n1, n2 = n2,
        estimate = NA_real_, statistic = NA_real_,
        p_value = NA_real_,
        method = method,
        stringsAsFactors = FALSE
      ))
    }
    est <- sm$coefficients[idx, "Estimate"]
    stat <- sm$coefficients[idx, "t value"]
    pv <- sm$coefficients[idx, "Pr(>|t|)"]
    data.frame(
      group1 = ref, group2 = target,
      n1 = n1, n2 = n2,
      estimate = unname(est),
      statistic = unname(stat),
      p_value = unname(pv),
      method = method,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, out)
  if (!is.data.frame(out) || nrow(out) == 0) {
    return(data.frame(
      group1 = character(0), group2 = character(0),
      n1 = integer(0), n2 = integer(0),
      estimate = numeric(0), statistic = numeric(0),
      p_value = numeric(0), p_adj = numeric(0),
      method = character(0),
      stringsAsFactors = FALSE
    ))
  }

  if (identical(p_adjust, "none")) {
    out$p_adj <- out$p_value
  } else {
    out$p_adj <- stats::p.adjust(out$p_value, method = p_adjust)
  }
  out
}

#' Dot Plot of Group Profiles
#'
#' @param summary_tbl Data frame with summary statistics (e.g. from
#'   \code{compare_groups()$summaries}).
#' @param se_col Character. Column name for standard error.
#' @param group_col Character. Column name for group labels.
#' @param var_col Character. Column name for variable/ROI names.
#' @param value_col Character. Column name for the mean (or plotted value).
#' @param plot_spec List. Optional plot customization (theme, sizes, colors).
#'
#' @return A \code{ggplot} object.
#' @export
plot_dot <- function(summary_tbl, se_col = "se", group_col = "groups",
                     var_col = "var", value_col = "mean",
                     plot_spec = list()) {
  .check_pkg("ggplot2", "plotting")
  plot_spec <- .as_plot_spec(plot_spec)
  point_size <- plot_spec$point_size %||% 2.2
  point_alpha <- plot_spec$point_alpha %||% 1
  point_shape <- plot_spec$point_shape %||% 16
  errorbar_width <- plot_spec$errorbar_width %||% 0.3
  x_angle <- plot_spec$x_text_angle %||% 45
  x_hjust <- plot_spec$x_text_hjust %||% 1

  df <- summary_tbl
  df$.var <- df[[var_col]]
  df$.value <- df[[value_col]]
  df$.group <- df[[group_col]]
  df$.ymin <- df[[value_col]] - df[[se_col]]
  df$.ymax <- df[[value_col]] + df[[se_col]]

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .var, y = .value, color = .group, group = .group)) +
    ggplot2::geom_point(size = point_size, alpha = point_alpha, shape = point_shape) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .ymin, ymax = .ymax), width = errorbar_width) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = x_angle, hjust = x_hjust)) +
    ggplot2::labs(x = NULL, y = "Value", color = group_col)
  p <- .apply_discrete_scales(p, plot_spec, use_color = TRUE, use_fill = FALSE)
  .apply_common_plot_spec(p, plot_spec)
}

#' Bar Plot of Group Profiles
#'
#' @param summary_tbl Data frame with summary statistics.
#' @param se_col Character. Column name for standard error.
#' @param group_col Character. Column name for group labels.
#' @param var_col Character. Column name for variable/ROI names.
#' @param value_col Character. Column name for the mean (or plotted value).
#' @param plot_spec List. Optional plot customization (theme, sizes, colors).
#'
#' @return A \code{ggplot} object.
#' @export
plot_bar <- function(summary_tbl, se_col = "se", group_col = "groups",
                     var_col = "var", value_col = "mean",
                     plot_spec = list()) {
  .check_pkg("ggplot2", "plotting")
  plot_spec <- .as_plot_spec(plot_spec)
  dodge_width <- plot_spec$dodge_width %||% 0.8
  bar_width <- plot_spec$bar_width %||% 0.8
  errorbar_width <- plot_spec$errorbar_width %||% 0.3

  df <- summary_tbl
  df$.var <- df[[var_col]]
  df$.value <- df[[value_col]]
  df$.group <- df[[group_col]]
  df$.ymin <- df[[value_col]] - df[[se_col]]
  df$.ymax <- df[[value_col]] + df[[se_col]]

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .var, y = .value, fill = .group)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = dodge_width), width = bar_width) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .ymin, ymax = .ymax),
                           position = ggplot2::position_dodge(width = dodge_width),
                           width = errorbar_width) +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = NULL, y = "Value", fill = group_col)
  p <- .apply_discrete_scales(p, plot_spec, use_color = FALSE, use_fill = TRUE)
  .apply_common_plot_spec(p, plot_spec)
}

#' Violin Plot (Raw Data, No Re-fitting)
#'
#' @param data Data frame.
#' @param x Character. X-axis variable (typically a grouping column).
#' @param y Character. Y-axis variable.
#' @param group_col Optional character. Column used for fill/group aesthetics.
#' @param plot_spec List. Optional plot customization (theme, sizes, colors).
#'
#' @return A \code{ggplot} object.
#' @export
plot_violin <- function(data, x, y, group_col = NULL, plot_spec = list()) {
  .check_pkg("ggplot2", "plotting")
  if (is.null(group_col)) group_col <- x
  plot_spec <- .as_plot_spec(plot_spec)
  violin_width <- plot_spec$violin_width %||% 0.6
  violin_alpha <- plot_spec$violin_alpha %||% 0.8
  box_width <- plot_spec$box_width %||% 0.08

  df <- data
  df$.x <- df[[x]]
  df$.y <- df[[y]]
  df$.group <- df[[group_col]]

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .x, y = .y, fill = .group)) +
    ggplot2::geom_violin(trim = TRUE, width = violin_width, alpha = violin_alpha) +
    ggplot2::geom_boxplot(width = box_width, outlier.shape = NA, alpha = plot_spec$box_alpha %||% 0.9) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = NULL, y = y, fill = group_col)
  p <- .apply_discrete_scales(p, plot_spec, use_color = FALSE, use_fill = TRUE)
  .apply_common_plot_spec(p, plot_spec)
}

#' Radar Plot for Group Profiles
#'
#' @param summary_tbl Data frame with summary statistics.
#' @param var_col Character. Column name for variable/ROI names.
#' @param value_col Character. Column name for the mean (or plotted value).
#' @param group_col Character. Column name for group labels.
#' @param vars Optional character vector to subset variables.
#' @param max_value Numeric. Max axis value used by radar chart.
#' @param min_value Numeric. Min axis value used by radar chart.
#' @param plot_spec List. Optional plot customization (theme, title, base_size).
#'
#' @return A named list of \code{ggplot} objects (one per group).
#' @export
plot_radar <- function(summary_tbl, var_col = "var", value_col = "mean",
                       group_col = "groups", vars = NULL,
                       max_value = 2, min_value = -2,
                       plot_spec = list()) {
  .check_pkg("fmsb", "radar plots")
  .check_pkg("ggplotify", "radar plots")
  plot_spec <- .as_plot_spec(plot_spec)
  if (!is.null(vars)) {
    summary_tbl <- summary_tbl[summary_tbl[[var_col]] %in% vars, , drop = FALSE]
  }
  groups <- unique(summary_tbl[[group_col]])
  plots <- list()
  for (g in groups) {
    df_g <- summary_tbl[summary_tbl[[group_col]] == g, , drop = FALSE]
    mat <- data.frame(
      max = max_value,
      min = min_value,
      setNames(as.list(df_g[[value_col]]), df_g[[var_col]])
    )
    mat <- as.data.frame(t(mat))
    p <- ggplotify::as.ggplot(function() fmsb::radarchart(mat))
    plots[[as.character(g)]] <- .apply_common_plot_spec(p, plot_spec)
  }
  plots
}

#' Quadrant Plot from Stat Table
#'
#' @param stat_tbl Data frame (typically a results table).
#' @param x_col Character. Column name for x-axis value.
#' @param y_col Character. Column name for y-axis value.
#' @param label_col Optional character. Column name used for point labels.
#' @param axis_max Numeric. Symmetric axis limit (±axis_max).
#' @param add_label Logical. Add repelled labels.
#' @param plot_spec List. Optional plot customization (theme, point size/alpha).
#'
#' @return A \code{ggplot} object.
#' @export
plot_quadrant <- function(stat_tbl, x_col, y_col, label_col = NULL,
                          axis_max = 1, add_label = FALSE,
                          plot_spec = list()) {
  .check_pkg("ggplot2", "plotting")
  plot_spec <- .as_plot_spec(plot_spec)
  point_size <- plot_spec$point_size %||% 2.5
  point_alpha <- plot_spec$point_alpha %||% 1

  df <- stat_tbl
  df$.x <- df[[x_col]]
  df$.y <- df[[y_col]]
  if (!is.null(label_col)) df$.label <- df[[label_col]]

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .x, y = .y)) +
    ggplot2::geom_point(size = point_size, alpha = point_alpha) +
    ggplot2::xlim(-axis_max, axis_max) +
    ggplot2::ylim(-axis_max, axis_max) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::theme_minimal()
  if (add_label && !is.null(label_col)) {
    .check_pkg("ggrepel", "labeling")
    p <- p + ggrepel::geom_label_repel(ggplot2::aes(label = .label))
  }
  .apply_common_plot_spec(p, plot_spec)
}

#' Brain Map Plot (ggseg-based)
#'
#' @description
#' Creates a brain map using \pkg{ggseg}. This function consumes a table (no
#' model fitting) and supports ggseg layout options via \code{plot_spec}.
#'
#' @param tbl Data frame with at least a label column (e.g. \code{"lh_bankssts"}).
#' @param atlas_spec List with fields: \code{atlas} (dk/desterieux/aseg),
#'   \code{label_col}, \code{value_col}.
#' @param plot_spec List. Common fields:
#'   \itemize{
#'     \item \code{layout}: \code{"dispersed"} (default) or \code{"stacked"}.
#'     \item \code{view}: optional view filter (e.g. \code{"lateral"}).
#'     \item \code{hemisphere}: optional hemi filter (e.g. \code{c("left","right")}).
#'     \item \code{p_col}, \code{p_max}: optional thresholding (set fill to NA if p > p_max).
#'     \item \code{limit}, \code{midpoint}, \code{legend_title}, \code{title}.
#'   }
#'
#' @export
plot_brain_map <- function(tbl, atlas_spec, plot_spec = list()) {
  .check_pkg("ggseg", "brain plotting")
  .check_pkg("ggplot2", "brain plotting")

  plot_spec <- .as_plot_spec(plot_spec)

  atlas <- atlas_spec$atlas %||% "dk"
  value_col <- atlas_spec$value_col %||% "value"
  label_col <- atlas_spec$label_col %||% "label"

  if (value_col %in% c("es", "es_value")) {
    if (!all(c("ref_level", "target_level", "es_method") %in% names(tbl))) {
      stop("Effect-size plotting requires ref_level, target_level, and es_method columns.", call. = FALSE)
    }
  }

  df <- tbl
  df$label <- df[[label_col]]

  # Optional thresholding at plot time
  p_col <- plot_spec$p_col %||% NULL
  p_max <- plot_spec$p_max %||% NULL
  if (!is.null(p_col) && !is.null(p_max)) {
    if (!p_col %in% names(df)) stop("p_col not found in tbl.", call. = FALSE)
    df[[value_col]][!(df[[p_col]] <= p_max)] <- NA_real_
  }

  # Use a stable column name for aesthetics to avoid deprecated aes_string().
  df$.value <- df[[value_col]]

  atlas_obj <- switch(
    atlas,
    dk = ggseg::dk,
    desterieux = {
      .check_pkg("ggsegDesterieux", "desterieux atlas")
      ggsegDesterieux::desterieux
    },
    aseg = ggseg::aseg,
    stop("Unknown atlas in atlas_spec.", call. = FALSE)
  )

  # Layout handling:
  # - position: "dispersed" (default) or "stacked"
  # - view: e.g., "lateral", "medial", "frontal" (passed to ggseg filter)
  position <- plot_spec$position %||% "dispersed"
  view <- plot_spec$view %||% NULL
  hemisphere <- plot_spec$hemisphere %||% NULL
  layout <- plot_spec$layout %||% NULL
  if (!is.null(layout)) {
    if (is.list(layout)) {
      if (!is.null(layout$position)) position <- layout$position
      if (!is.null(layout$view)) view <- layout$view
      if (!is.null(layout$hemisphere)) hemisphere <- layout$hemisphere
    } else if (is.character(layout) && length(layout) == 1) {
      if (layout %in% c("stacked", "dispersed")) {
        position <- layout
      } else if (layout == "default") {
        position <- "dispersed"
        view <- NULL
      } else if (is.null(view)) {
        view <- layout
      }
    }
  }

  limit <- plot_spec$limit %||% NULL
  midpoint <- plot_spec$midpoint %||% 0
  color_low <- plot_spec$color_low %||% "steelblue1"
  color_high <- plot_spec$color_high %||% "firebrick1"
  color_mid <- plot_spec$color_mid %||% "white"
  na_color <- plot_spec$na_color %||% "grey90"
  legend_title <- plot_spec$legend_title %||% "Stat"
  title <- plot_spec$title %||% ""
  hide_legend <- isTRUE(plot_spec$hide_legend)

  if (is.null(limit)) {
    mx <- max(abs(df[[value_col]]), na.rm = TRUE)
    if (!is.finite(mx) || mx == 0) mx <- 1
    limit <- c(-mx, mx)
  }

  p <- ggseg::ggseg(
    df,
    atlas = atlas_obj,
    position = position,
    view = view,
    hemisphere = hemisphere,
    mapping = ggplot2::aes(fill = .value),
    colour = plot_spec$border_color %||% "grey40",
    size = plot_spec$border_size %||% 0.1
  ) +
    ggplot2::scale_fill_gradient2(
      midpoint = midpoint,
      low = color_low, mid = color_mid, high = color_high,
      limits = limit,
      na.value = na_color
    ) +
    ggplot2::labs(fill = legend_title, title = title)

  # Use ggseg-native themes only for brain maps. Avoid ggplot2 themes like
  # theme_minimal() here because they can re-introduce panel grid lines that
  # look like vertical separators between ggseg facets.
  ggseg_theme <- plot_spec$ggseg_theme %||% "brain2"
  p <- p + switch(
    ggseg_theme,
    brain = ggseg::theme_brain(),
    brain2 = ggseg::theme_brain2(),
    darkbrain = ggseg::theme_darkbrain(),
    custombrain = ggseg::theme_custombrain(),
    stop("Unknown plot_spec$ggseg_theme (use: brain/brain2/darkbrain/custombrain).", call. = FALSE)
  )
  # Extra hardening: blank any remaining grids to prevent separator artifacts.
  p <- p + ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid = ggplot2::element_blank()
  )

  if (hide_legend) p <- p + ggplot2::theme(legend.position = "none")
  if (!is.null(plot_spec$legend_position)) {
    p <- p + ggplot2::theme(legend.position = plot_spec$legend_position)
  }
  if (!is.null(plot_spec$base_size)) {
    p <- p + ggplot2::theme(text = ggplot2::element_text(size = plot_spec$base_size))
  }
  if (!is.null(plot_spec$subtitle)) p <- p + ggplot2::labs(subtitle = plot_spec$subtitle)
  if (!is.null(plot_spec$caption)) p <- p + ggplot2::labs(caption = plot_spec$caption)
  p
}

#' Brain Map from a Results Table (ROI Names -> ggseg Labels)
#'
#' @description
#' Convenience wrapper around \code{\link{plot_brain_map}} that converts
#' FreeSurfer-style ROI names (e.g. \code{"L_superiorfrontal_thickavg"}) into
#' ggseg labels (e.g. \code{"lh_superiorfrontal"}) and then renders a brain map.
#'
#' @param res_tbl Data frame (typically a results table).
#' @param roi_col Character. Column containing ROI names. Default: \code{"var"}.
#' @param value_col Character. Column containing the value to plot.
#' @param atlas Character. One of \code{"dk"}, \code{"desterieux"}, \code{"aseg"}.
#' @param roi_drop_suffix Optional suffix removed from ROI names before mapping.
#' @param plot_spec List. Passed through to \code{\link{plot_brain_map}}.
#'
#' @return A \code{ggplot} object.
#' @export
plot_brain_map_results <- function(res_tbl,
                                   roi_col = "var",
                                   value_col = "es_value",
                                   atlas = c("dk", "desterieux", "aseg"),
                                   roi_drop_suffix = "_thickavg",
                                   plot_spec = list()) {
  atlas <- match.arg(atlas)
  .check_pkg("ggplot2", "plotting")
  stopifnot(is.data.frame(res_tbl))
  if (!roi_col %in% names(res_tbl)) stop("roi_col not found in res_tbl.", call. = FALSE)
  if (!value_col %in% names(res_tbl)) stop("value_col not found in res_tbl.", call. = FALSE)

  df <- res_tbl
  df$label <- .roi_to_ggseg_label(df[[roi_col]], drop_suffix = roi_drop_suffix)

  plot_brain_map(
    df,
    atlas_spec = list(atlas = atlas, label_col = "label", value_col = value_col),
    plot_spec = plot_spec
  )
}

#' QC Plot: P-Value Histogram
#'
#' @param res_tbl Data frame (results table).
#' @param p_col Character. Column name containing p-values (e.g. \code{"p_adj"}).
#' @param nbins Integer. Histogram bin count.
#' @param title Character. Plot title.
#' @param plot_spec List. Optional plot customization (theme, fill/alpha).
#'
#' @return A \code{ggplot} object.
#' @export
plot_qc_pvalue_hist <- function(res_tbl, p_col = "p_adj", nbins = 30,
                                title = "QC: p-value distribution",
                                plot_spec = list()) {
  .check_pkg("ggplot2", "plotting")
  stopifnot(is.data.frame(res_tbl))
  if (!p_col %in% names(res_tbl)) stop("p_col not found in res_tbl.", call. = FALSE)
  plot_spec <- .as_plot_spec(plot_spec)
  fill <- plot_spec$fill %||% "grey"
  alpha <- plot_spec$alpha %||% 0.85

  df <- data.frame(p = res_tbl[[p_col]])
  p <- ggplot2::ggplot(df, ggplot2::aes(x = p)) +
    ggplot2::geom_histogram(bins = nbins, fill = fill, alpha = alpha) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = p_col, y = "Count", title = title)
  .apply_common_plot_spec(p, plot_spec)
}

#' Top ROIs Bar Plot
#'
#' @param res_tbl Data frame (results table).
#' @param metric Character. Column used for ranking/plotting (absolute value).
#' @param n Integer. Number of rows to plot.
#' @param label_col Character. Column used as labels.
#' @param title Character. Plot title.
#' @param plot_spec List. Optional plot customization (theme, gradient colors).
#'
#' @return A \code{ggplot} object.
#' @export
plot_top_rois <- function(res_tbl, metric = "es_value", n = 20,
                          label_col = "var", title = "Top ROIs",
                          plot_spec = list()) {
  .check_pkg("ggplot2", "plotting")
  stopifnot(is.data.frame(res_tbl))
  if (!metric %in% names(res_tbl)) stop("metric not found in res_tbl.", call. = FALSE)
  if (!label_col %in% names(res_tbl)) stop("label_col not found in res_tbl.", call. = FALSE)
  n <- as.integer(n)
  if (!is.finite(n) || n < 1) stop("n must be >= 1.", call. = FALSE)

  ord <- order(abs(res_tbl[[metric]]), decreasing = TRUE)
  df <- res_tbl[ord, , drop = FALSE]
  df <- df[seq_len(min(n, nrow(df))), , drop = FALSE]
  df$label <- df[[label_col]]
  df$label <- factor(df$label, levels = rev(df$label))

  plot_spec <- .as_plot_spec(plot_spec)
  low <- plot_spec$color_low %||% "steelblue1"
  mid <- plot_spec$color_mid %||% "white"
  high <- plot_spec$color_high %||% "firebrick1"

  df$.metric <- df[[metric]]
  p <- ggplot2::ggplot(df, ggplot2::aes(x = label, y = .metric, fill = .metric)) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_gradient2(low = low, mid = mid, high = high, midpoint = 0) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = NULL, y = metric, title = title, fill = metric)
  .apply_common_plot_spec(p, plot_spec)
}

#' AUC Histogram Plot
#'
#' @param auc_null Numeric vector. Null AUC distribution.
#' @param auc_true Optional numeric. Observed/true AUC.
#' @param nbins Integer. Histogram bin count.
#' @param xlim Optional numeric vector of length 2. X limits.
#' @param ylim Optional numeric vector of length 2. Y limits.
#' @param line_color Character. Color for the observed AUC line.
#' @param plot_spec List. Optional plot customization (theme, fill/alpha).
#'
#' @return A \code{ggplot} object.
#' @export
plot_auc_hist <- function(auc_null, auc_true = NULL, nbins = 30,
                          xlim = NULL, ylim = NULL, line_color = "steelblue1",
                          plot_spec = list()) {
  .check_pkg("ggplot2", "plotting")
  plot_spec <- .as_plot_spec(plot_spec)
  df <- data.frame(auc = auc_null)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = auc)) +
    ggplot2::geom_histogram(bins = nbins, fill = plot_spec$fill %||% "grey",
                            alpha = plot_spec$alpha %||% 0.8) +
    ggplot2::geom_vline(xintercept = mean(auc_null, na.rm = TRUE),
                        color = "grey3", linetype = "longdash")
  if (!is.null(auc_true)) {
    p <- p + ggplot2::geom_vline(xintercept = auc_true, color = line_color)
  }
  if (!is.null(xlim)) p <- p + ggplot2::xlim(xlim)
  if (!is.null(ylim)) p <- p + ggplot2::ylim(ylim)
  p <- p + ggplot2::theme_minimal() + ggplot2::labs(x = "AUC", y = "Count")
  .apply_common_plot_spec(p, plot_spec)
}

#' Raincloud Plot for a Selected ROI (Half-Violin + Points)
#'
#' @description
#' Raincloud-style distribution plot: a half-violin density ("cloud") plus
#' jittered points ("rain"). Supports 2+ groups and optional pairwise
#' significance markers.
#'
#' \code{value = "marginal"} produces covariate-adjusted ROI values derived from
#' a linear model \code{ROI ~ group + covariates}, where covariate terms are
#' removed and centered at their average contribution (group effects retained).
#'
#' @param data Data frame.
#' @param roi Character. ROI column name.
#' @param group_col Character. Grouping column (>= 2 levels).
#' @param covariates Optional character vector of covariate column names.
#' @param value Character. \code{"raw"} or \code{"marginal"}.
#' @param orientation Character. \code{"vertical"} (default) or \code{"horizontal"}.
#' @param style Character. Currently only \code{"half_violin"} is supported.
#' @param pairwise Optional list enabling pairwise comparisons + markers. Fields:
#'   \itemize{
#'     \item \code{method}: \code{"t_test"}, \code{"wilcox"}, or \code{"lm"}.
#'     \item \code{p_adjust}: \code{"holm"}, \code{"fdr"}, \code{"bonferroni"}, \code{"none"}.
#'     \item \code{pairs}: optional list of pairs (each length-2 character vector).
#'     \item \code{label}: \code{"stars"} (default) or \code{"p_adj"}.
#'     \item \code{p_col}: \code{"p_adj"} (default) or \code{"p_value"}.
#'     \item \code{show_ns}: show non-significant labels (default: FALSE).
#'   }
#' @param plot_spec List. Customization options (theme, sizes, colors, etc.).
#'
#' @return A \code{ggplot} object.
#' @export
plot_raincloud_roi <- function(data, roi, group_col,
                               covariates = NULL,
                               value = c("raw", "marginal"),
                               orientation = c("vertical", "horizontal"),
                               style = c("half_violin"),
                               pairwise = NULL,
                               plot_spec = list()) {
  value <- match.arg(value)
  orientation <- match.arg(orientation)
  style <- match.arg(style)
  plot_spec <- .as_plot_spec(plot_spec)

  .check_pkg("ggplot2", "plotting")
  .check_pkg("gghalves", "raincloud plots")
  stopifnot(is.data.frame(data))
  .ensure_cols(data, c(roi, group_col, covariates), context = "plot_raincloud_roi")

  # Only require complete cases for what we actually plot.
  cols_plot <- c(roi, group_col)
  if (identical(value, "marginal") && !is.null(covariates) && length(covariates) > 0) {
    cols_plot <- c(cols_plot, covariates)
  }
  keep <- stats::complete.cases(data[, cols_plot, drop = FALSE])
  d <- data[keep, , drop = FALSE]
  if (nrow(d) == 0) stop("No complete cases for raincloud plot.", call. = FALSE)

  # Group ordering: prefer explicit plot_spec$group_levels, else preserve factor levels.
  g0 <- d[[group_col]]
  if (!is.null(plot_spec$group_levels)) {
    levs <- as.character(plot_spec$group_levels)
    if (any(!levs %in% unique(as.character(g0)))) {
      stop("plot_spec$group_levels contains levels not found in data.", call. = FALSE)
    }
    d[[group_col]] <- factor(as.character(g0), levels = levs)
  } else if (is.factor(g0)) {
    d[[group_col]] <- droplevels(g0)
  } else {
    levs <- unique(as.character(g0))
    d[[group_col]] <- factor(as.character(g0), levels = levs)
  }

  if (nlevels(d[[group_col]]) < 2) {
    stop("plot_raincloud_roi requires >= 2 group levels.", call. = FALSE)
  }

  y_plot <- d[[roi]]
  if (identical(value, "marginal") && !is.null(covariates) && length(covariates) > 0) {
    rhs <- paste(c(.quote_var(group_col), vapply(covariates, .quote_var, character(1))), collapse = " + ")
    f <- stats::as.formula(paste(.quote_var(roi), "~", rhs))
    fit <- stats::lm(f, data = d)
    y_plot <- .adjust_outcome_keep_terms(fit, d, y_col = roi, keep_terms = group_col, center = TRUE)
  }

  df_plot <- data.frame(
    group = d[[group_col]],
    value = y_plot,
    stringsAsFactors = FALSE
  )

  # Raincloud geometry knobs
  violin_alpha <- plot_spec$violin_alpha %||% 0.75
  violin_width <- plot_spec$violin_width %||% 0.8
  violin_scale <- plot_spec$violin_scale %||% "width"
  violin_side <- plot_spec$violin_side %||% "l"
  violin_nudge <- plot_spec$violin_nudge %||% 0.05

  point_size <- plot_spec$point_size %||% 1.6
  point_alpha <- plot_spec$point_alpha %||% 0.55
  point_side <- plot_spec$point_side %||% "r"
  jitter_width <- plot_spec$jitter_width %||% 0.1
  range_scale <- plot_spec$point_range_scale %||% 0.7

  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = group, y = value, fill = group, color = group)) +
    gghalves::geom_half_violin(
      side = violin_side,
      nudge = violin_nudge,
      width = violin_width,
      alpha = violin_alpha,
      trim = TRUE,
      scale = violin_scale
    ) +
    gghalves::geom_half_point(
      side = point_side,
      size = point_size,
      alpha = point_alpha,
      transformation = ggplot2::position_jitter(width = jitter_width, height = 0),
      range_scale = range_scale
    ) +
    ggplot2::labs(
      x = NULL,
      y = if (identical(value, "marginal")) paste0(roi, " (marginal)") else roi,
      title = plot_spec$title %||% sprintf("Raincloud: %s by %s (%s)", roi, group_col, value)
    )

  # Optional central tendency marker (mean +/- normal CI)
  if (isTRUE(plot_spec$show_summary %||% TRUE)) {
    p <- p + ggplot2::stat_summary(
      fun.data = ggplot2::mean_cl_normal,
      geom = "pointrange",
      linewidth = plot_spec$summary_linewidth %||% 0.5,
      size = plot_spec$summary_size %||% 0.2,
      colour = plot_spec$summary_color %||% "black"
    )
  }

  # Pairwise comparisons + markers
  if (!is.null(pairwise)) {
    pw <- if (isTRUE(pairwise)) list() else pairwise
    if (!is.list(pw)) stop("pairwise must be NULL, TRUE, or a list.", call. = FALSE)
    pw_method <- pw$method %||% "lm"
    pw_p_adjust <- pw$p_adjust %||% "holm"
    pw_pairs <- pw$pairs %||% NULL
    pw_label <- pw$label %||% "stars"
    pw_p_col <- pw$p_col %||% "p_adj"

    tests <- pairwise_tests_roi(
      data = data,
      roi = roi,
      group_col = group_col,
      covariates = covariates,
      value = value,
      method = pw_method,
      p_adjust = pw_p_adjust,
      pairs = pw_pairs
    )

    if (nrow(tests) > 0 && pw_p_col %in% names(tests) && any(is.finite(tests[[pw_p_col]]))) {
      p_vals <- tests[[pw_p_col]]
      labels <- if (identical(pw_label, "p_adj")) {
        sprintf("p=%s", format(p_vals, digits = 2, scientific = TRUE))
      } else {
        vapply(p_vals, .p_to_stars, character(1))
      }

      show_ns <- isTRUE(pw$show_ns %||% FALSE)
      keep_idx <- which(!is.na(labels) & (show_ns | labels != "ns"))
      tests <- tests[keep_idx, , drop = FALSE]
      labels <- labels[keep_idx]

      if (nrow(tests) > 0) {
        y_max <- max(df_plot$value, na.rm = TRUE)
        rng <- diff(range(df_plot$value, na.rm = TRUE))
        if (!is.finite(rng) || rng == 0) rng <- 1
        step <- plot_spec$pairwise_step %||% (0.06 * rng)
        tick <- plot_spec$pairwise_tick %||% (0.02 * rng)
        text_off <- plot_spec$pairwise_text_offset %||% (0.02 * rng)

        levs <- levels(df_plot$group)
        x1 <- match(as.character(tests$group1), levs)
        x2 <- match(as.character(tests$group2), levs)
        y <- y_max + seq_len(nrow(tests)) * step

        p <- p + ggplot2::coord_cartesian(clip = "off") +
          ggplot2::expand_limits(y = max(y + text_off, na.rm = TRUE)) +
          ggplot2::annotate("segment", x = x1, xend = x2, y = y, yend = y,
                            linewidth = plot_spec$pairwise_linewidth %||% 0.5) +
          ggplot2::annotate("segment", x = x1, xend = x1, y = y, yend = y - tick,
                            linewidth = plot_spec$pairwise_linewidth %||% 0.5) +
          ggplot2::annotate("segment", x = x2, xend = x2, y = y, yend = y - tick,
                            linewidth = plot_spec$pairwise_linewidth %||% 0.5) +
          ggplot2::annotate("text", x = (x1 + x2) / 2, y = y + text_off, label = labels,
                            size = plot_spec$pairwise_text_size %||% 3)
      }
    }
  }

  p <- .apply_discrete_scales(p, plot_spec, use_color = TRUE, use_fill = TRUE)
  p <- .apply_common_plot_spec(p, plot_spec)
  if (identical(orientation, "horizontal")) p <- p + ggplot2::coord_flip()
  p
}

#' Scatter Plot for Association (Selected ROI)
#'
#' @description
#' Scatter plot of a selected ROI versus a predictor variable (e.g. Age),
#' optionally adjusting the ROI values for covariates.
#'
#' \code{value = "marginal"} produces covariate-adjusted ROI values derived from
#' a linear model \code{ROI ~ x + covariates}, where covariate terms (excluding
#' \code{x}) are removed and centered at their average.
#'
#' @param data Data frame.
#' @param roi Character. ROI column name.
#' @param x Character. Predictor variable (x-axis).
#' @param covariates Optional character vector of covariates to adjust for.
#' @param color Optional character. Column used for point colors (e.g. Sex).
#' @param value Character. \code{"raw"} or \code{"marginal"}.
#' @param plot_spec List. Customization options (theme, sizes, colors).
#'
#' @return A \code{ggplot} object.
#' @export
plot_scatter_assoc <- function(data, roi, x,
                               covariates = NULL,
                               color = NULL,
                               value = c("raw", "marginal"),
                               plot_spec = list()) {
  value <- match.arg(value)
  plot_spec <- .as_plot_spec(plot_spec)

  .check_pkg("ggplot2", "plotting")
  stopifnot(is.data.frame(data))

  cols_needed <- unique(c(roi, x, covariates, color))
  cols_needed <- cols_needed[!is.na(cols_needed) & nzchar(cols_needed)]
  .ensure_cols(data, cols_needed, context = "plot_scatter_assoc")

  keep <- stats::complete.cases(data[, cols_needed, drop = FALSE])
  d <- data[keep, , drop = FALSE]
  if (nrow(d) == 0) stop("No complete cases for scatter plot.", call. = FALSE)

  y_raw <- d[[roi]]
  y_plot <- y_raw

  rhs <- paste(c(.quote_var(x), vapply(covariates %||% character(0), .quote_var, character(1))), collapse = " + ")
  f <- stats::as.formula(paste(.quote_var(roi), "~", rhs))

  if (identical(value, "marginal") && !is.null(covariates) && length(covariates) > 0) {
    fit <- stats::lm(f, data = d)
    y_plot <- .adjust_outcome_keep_terms(fit, d, y_col = roi, keep_terms = x, center = TRUE)
  }

  df_plot <- data.frame(
    x = d[[x]],
    y = y_plot,
    stringsAsFactors = FALSE
  )
  if (!is.null(color)) df_plot$color <- d[[color]]

  point_size <- plot_spec$point_size %||% 1.6
  point_alpha <- plot_spec$point_alpha %||% 0.6
  se <- plot_spec$se %||% TRUE

  mapping <- if (is.null(color)) {
    ggplot2::aes(x = x, y = y)
  } else {
    ggplot2::aes(x = x, y = y, color = color)
  }

  p <- ggplot2::ggplot(df_plot, mapping) +
    ggplot2::geom_point(size = point_size, alpha = point_alpha) +
    ggplot2::geom_smooth(method = "lm", se = se, linewidth = plot_spec$line_size %||% 0.8) +
    ggplot2::labs(
      x = x,
      y = if (identical(value, "marginal")) paste0(roi, " (marginal)") else roi,
      title = plot_spec$title %||% sprintf("Association: %s vs %s (%s)", roi, x, value)
    )

  if (!is.null(color) && !is.null(plot_spec$colors)) {
    p <- p + ggplot2::scale_color_manual(values = plot_spec$colors)
  }
  .apply_common_plot_spec(p, plot_spec)
}
