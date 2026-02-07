# PET correlation plotting utilities.

.pet_plot_theme <- function(theme = c("minimal", "classic", "bw"), base_size = 11) {
  theme <- match.arg(theme)
  switch(
    theme,
    minimal = ggplot2::theme_minimal(base_size = base_size),
    classic = ggplot2::theme_classic(base_size = base_size),
    bw = ggplot2::theme_bw(base_size = base_size)
  )
}

.pet_plot_result_check <- function(pet_corr_obj) {
  if (!inherits(pet_corr_obj, "pet_corr_result")) {
    .pet_stop("E_PET_PLOT_INPUT", "`pet_corr_obj` must be a `pet_corr_result`.")
  }
  if (is.null(pet_corr_obj$results) || !is.data.frame(pet_corr_obj$results)) {
    .pet_stop("E_PET_PLOT_INPUT", "`pet_corr_obj$results` is missing or invalid.")
  }
  if (is.null(pet_corr_obj$inputs$brain) || is.null(pet_corr_obj$inputs$pet)) {
    .pet_stop("E_PET_PLOT_INPUT", "`pet_corr_obj` does not include `inputs$brain` and `inputs$pet`.")
  }
  TRUE
}

.pet_plot_annotation <- function(position = c("top_left", "top_right", "bottom_left", "bottom_right")) {
  position <- match.arg(position)
  switch(
    position,
    top_left = list(x = -Inf, y = Inf, hjust = -0.05, vjust = 1.15),
    top_right = list(x = Inf, y = Inf, hjust = 1.05, vjust = 1.15),
    bottom_left = list(x = -Inf, y = -Inf, hjust = -0.05, vjust = -0.1),
    bottom_right = list(x = Inf, y = -Inf, hjust = 1.05, vjust = -0.1)
  )
}

.pet_plot_filter_results <- function(pet_corr_obj, brain_map = NULL) {
  res <- pet_corr_obj$results
  if (!is.null(brain_map)) {
    brain_map <- as.character(brain_map)
    bad <- setdiff(brain_map, unique(res$brain_map))
    if (length(bad) > 0) {
      .pet_stop("E_PET_PLOT_MAP", sprintf("Unknown `brain_map`: %s", paste(bad, collapse = ", ")))
    }
    res <- res[res$brain_map %in% brain_map, , drop = FALSE]
  }
  if (nrow(res) == 0) {
    .pet_stop("E_PET_PLOT_EMPTY", "No rows available for plotting after filtering.")
  }
  res
}

#' Scatter Plot for One Brain-PET Pair from `pet_corr()`
#'
#' @description
#' Plots ROI-level scatter values for a selected brain map and PET map from a
#' \code{pet_corr_result}, with optional regression line and mandatory
#' annotation of \code{r} and \code{p_spin}.
#'
#' @param pet_corr_obj A \code{pet_corr_result} object.
#' @param brain_map Character. Brain map column name to plot.
#' @param pet_map Character. PET map column name to plot.
#' @param add_smooth Logical. Add regression line.
#' @param smooth_method Character method passed to \code{geom_smooth()}.
#' @param smooth_se Logical. Show uncertainty band for smooth line.
#' @param point_color Point color.
#' @param point_size Point size.
#' @param point_alpha Point alpha.
#' @param line_color Smoothing line color.
#' @param line_size Smoothing line size.
#' @param annotation_position Corner for annotation text.
#' @param annotation_digits Digits for \code{r} formatting.
#' @param annotation_size Annotation label size.
#' @param annotation_fill Annotation label fill.
#' @param annotation_alpha Annotation label alpha.
#' @param title Optional title.
#' @param subtitle Optional subtitle.
#' @param xlab Optional x-axis label.
#' @param ylab Optional y-axis label.
#' @param theme Plot theme: \code{"minimal"}, \code{"classic"}, \code{"bw"}.
#' @param base_size Base font size.
#'
#' @return A \code{ggplot} object.
#' @export
#'
#' @examples
#' data(pet_maps_dk68)
#' perm_id <- replicate(20, sample.int(nrow(pet_maps_dk68)), simplify = "matrix")
#' res <- pet_corr(
#'   brain = pet_maps_dk68[, c("D1", "DAT")],
#'   pet_select = c("D1", "DAT"),
#'   perm_id = perm_id,
#'   n_perm = 20
#' )
#' p <- pet_plot_scatter(res, brain_map = "D1", pet_map = "D1")
pet_plot_scatter <- function(pet_corr_obj,
                             brain_map,
                             pet_map,
                             add_smooth = TRUE,
                             smooth_method = "lm",
                             smooth_se = FALSE,
                             point_color = "#2C7FB8",
                             point_size = 2,
                             point_alpha = 0.8,
                             line_color = "#D95F0E",
                             line_size = 0.8,
                             annotation_position = c("top_left", "top_right", "bottom_left", "bottom_right"),
                             annotation_digits = 3,
                             annotation_size = 3.6,
                             annotation_fill = "white",
                             annotation_alpha = 0.9,
                             title = NULL,
                             subtitle = NULL,
                             xlab = NULL,
                             ylab = NULL,
                             theme = c("minimal", "classic", "bw"),
                             base_size = 11) {
  .check_pkg("ggplot2", "PET plotting")
  .pet_plot_result_check(pet_corr_obj)
  theme <- match.arg(theme)

  brain_map <- as.character(brain_map)
  pet_map <- as.character(pet_map)

  bmat <- pet_corr_obj$inputs$brain
  pmat <- pet_corr_obj$inputs$pet
  if (!brain_map %in% colnames(bmat)) {
    .pet_stop("E_PET_PLOT_MAP", sprintf("`brain_map` not found: %s", brain_map))
  }
  if (!pet_map %in% colnames(pmat)) {
    .pet_stop("E_PET_PLOT_MAP", sprintf("`pet_map` not found: %s", pet_map))
  }

  df <- data.frame(
    brain = bmat[, brain_map],
    pet = pmat[, pet_map],
    stringsAsFactors = FALSE
  )
  keep <- stats::complete.cases(df)
  df <- df[keep, , drop = FALSE]
  if (nrow(df) == 0) {
    .pet_stop("E_PET_PLOT_EMPTY", "No complete ROI pairs available for scatter plot.")
  }

  stat_row <- pet_corr_obj$results[
    pet_corr_obj$results$brain_map == brain_map & pet_corr_obj$results$pet_map == pet_map,
    ,
    drop = FALSE
  ]
  if (nrow(stat_row) == 0) {
    .pet_stop("E_PET_PLOT_MAP", "Requested brain/PET pair not found in `pet_corr_obj$results`.")
  }
  stat_row <- stat_row[1, , drop = FALSE]

  annotation_position <- .pet_plot_annotation(annotation_position)
  label <- sprintf(
    "r = %.*f\np_spin = %.3g",
    as.integer(annotation_digits),
    stat_row$r_obs,
    stat_row$p_spin
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = brain, y = pet)) +
    ggplot2::geom_point(
      color = point_color,
      size = point_size,
      alpha = point_alpha
    ) +
    ggplot2::annotate(
      "label",
      x = annotation_position$x,
      y = annotation_position$y,
      hjust = annotation_position$hjust,
      vjust = annotation_position$vjust,
      label = label,
      size = annotation_size,
      fill = annotation_fill,
      alpha = annotation_alpha,
      label.size = 0.2
    ) +
    ggplot2::labs(
      x = xlab %||% brain_map,
      y = ylab %||% pet_map,
      title = title %||% sprintf("PET correlation scatter: %s vs %s", brain_map, pet_map),
      subtitle = subtitle
    ) +
    .pet_plot_theme(theme = theme, base_size = base_size)

  if (isTRUE(add_smooth)) {
    p <- p + ggplot2::geom_smooth(
      method = smooth_method,
      se = smooth_se,
      color = line_color,
      linewidth = line_size
    )
  }

  p
}

#' Bar/Dumbbell Plot for PET Correlation Results
#'
#' @description
#' Visualizes PET correlation outputs from \code{pet_corr()} across one or more
#' brain maps. Supports bar or dumbbell style.
#'
#' @param pet_corr_obj A \code{pet_corr_result} object.
#' @param brain_map Optional character vector to filter brain maps.
#' @param metric Which metric to plot on value axis:
#' \code{"r_obs"}, \code{"p_spin"}, or \code{"p_fdr"}.
#' @param style Plot style: \code{"bar"} or \code{"dumbbell"}.
#' @param order_by Ordering rule for PET maps.
#' @param top_n Optional integer. Keep top N PET maps per brain map by
#' \code{|r_obs|}.
#' @param sig_threshold Threshold used for alpha highlighting.
#' @param fill_low Low color for diverging scale.
#' @param fill_mid Mid color for diverging scale.
#' @param fill_high High color for diverging scale.
#' @param point_color Base point/segment color in dumbbell mode.
#' @param alpha_sig Alpha for significant rows.
#' @param alpha_nonsig Alpha for non-significant rows.
#' @param title Optional title.
#' @param xlab Optional x-axis label.
#' @param ylab Optional y-axis label.
#' @param legend_position Legend position.
#' @param facet_scales Facet scale behavior for multi-map plots.
#' @param theme Plot theme.
#' @param base_size Base font size.
#'
#' @return A \code{ggplot} object.
#' @export
#'
#' @examples
#' data(pet_maps_dk68)
#' perm_id <- replicate(20, sample.int(nrow(pet_maps_dk68)), simplify = "matrix")
#' res <- pet_corr(
#'   brain = pet_maps_dk68[, c("D1", "DAT")],
#'   pet_select = c("D1", "D2", "DAT", "MOR"),
#'   perm_id = perm_id,
#'   n_perm = 20
#' )
#' p <- pet_plot_bar(res, style = "bar")
pet_plot_bar <- function(pet_corr_obj,
                         brain_map = NULL,
                         metric = c("r_obs", "p_spin", "p_fdr"),
                         style = c("bar", "dumbbell"),
                         order_by = c("abs_r", "r_obs", "p_spin"),
                         top_n = NULL,
                         sig_threshold = 0.05,
                         fill_low = "#4575B4",
                         fill_mid = "white",
                         fill_high = "#D73027",
                         point_color = "#1F2937",
                         alpha_sig = 1,
                         alpha_nonsig = 0.45,
                         title = NULL,
                         xlab = NULL,
                         ylab = NULL,
                         legend_position = "right",
                         facet_scales = "fixed",
                         theme = c("minimal", "classic", "bw"),
                         base_size = 11) {
  .check_pkg("ggplot2", "PET plotting")
  .pet_plot_result_check(pet_corr_obj)
  metric <- match.arg(metric)
  style <- match.arg(style)
  order_by <- match.arg(order_by)
  theme <- match.arg(theme)

  df <- .pet_plot_filter_results(pet_corr_obj, brain_map = brain_map)

  if (!is.null(top_n)) {
    top_n <- as.integer(top_n)
    if (!is.finite(top_n) || top_n < 1) {
      .pet_stop("E_PET_PLOT_TOPN", "`top_n` must be a positive integer.")
    }
    keep_idx <- unlist(
      lapply(split(seq_len(nrow(df)), df$brain_map), function(idx) {
        ord <- order(abs(df$r_obs[idx]), decreasing = TRUE)
        idx[ord[seq_len(min(top_n, length(idx)))]]
      }),
      use.names = FALSE
    )
    df <- df[keep_idx, , drop = FALSE]
  }
  if (nrow(df) == 0) {
    .pet_stop("E_PET_PLOT_EMPTY", "No rows available for `pet_plot_bar()` after filtering.")
  }

  pet_score <- switch(
    order_by,
    abs_r = tapply(abs(df$r_obs), df$pet_map, max, na.rm = TRUE),
    r_obs = tapply(df$r_obs, df$pet_map, mean, na.rm = TRUE),
    p_spin = -tapply(df$p_spin, df$pet_map, mean, na.rm = TRUE)
  )
  pet_levels <- names(sort(pet_score, decreasing = TRUE))
  df$pet_map <- factor(df$pet_map, levels = rev(pet_levels))
  df$.metric <- df[[metric]]
  df$.sig <- ifelse(df$p_spin <= sig_threshold, "sig", "nonsig")

  single_brain <- length(unique(df$brain_map)) == 1
  if (style == "bar") {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = pet_map, y = .metric, alpha = .sig))
    if (metric == "r_obs") {
      p <- p + ggplot2::geom_col(ggplot2::aes(fill = .metric), width = 0.75)
      p <- p + ggplot2::scale_fill_gradient2(
        low = fill_low, mid = fill_mid, high = fill_high, midpoint = 0
      )
    } else {
      p <- p + ggplot2::geom_col(fill = fill_high, width = 0.75)
    }
    p <- p + ggplot2::coord_flip()
  } else {
    p <- ggplot2::ggplot(df, ggplot2::aes(y = pet_map, x = .metric, alpha = .sig))
    p <- p +
      ggplot2::geom_segment(
        ggplot2::aes(x = 0, xend = .metric, yend = pet_map),
        color = "grey70",
        linewidth = 0.6
      )
    if (metric == "r_obs") {
      p <- p + ggplot2::geom_point(ggplot2::aes(color = .metric), size = 2.5)
      p <- p + ggplot2::scale_color_gradient2(
        low = fill_low, mid = fill_mid, high = fill_high, midpoint = 0
      )
    } else {
      p <- p + ggplot2::geom_point(color = point_color, size = 2.5)
    }
  }

  p <- p +
    ggplot2::scale_alpha_manual(values = c(sig = alpha_sig, nonsig = alpha_nonsig), name = NULL) +
    ggplot2::labs(
      x = xlab %||% if (style == "bar") "PET map" else metric,
      y = ylab %||% if (style == "bar") metric else "PET map",
      title = title %||% sprintf("PET correlation %s plot (%s)", style, metric),
      fill = metric,
      color = metric
    ) +
    .pet_plot_theme(theme = theme, base_size = base_size) +
    ggplot2::theme(legend.position = legend_position)

  if (!single_brain) {
    p <- p + ggplot2::facet_wrap(~brain_map, scales = facet_scales)
  }
  p
}

#' Radar Plot for Multi-PET Correlation Profiles
#'
#' @description
#' Plots PET-correlation profiles from \code{pet_corr()} in radar format for one
#' or more brain maps.
#'
#' @param pet_corr_obj A \code{pet_corr_result} object.
#' @param brain_map Optional character vector to filter brain maps.
#' @param metric Value to plot: \code{"r_obs"} or \code{"minus_log10_p_spin"}.
#' @param pet_order Optional character vector defining PET map order.
#' @param line_size Polygon border width.
#' @param point_size Point size.
#' @param fill_alpha Polygon fill alpha.
#' @param palette Optional named vector of colors by brain map.
#' @param radial_limits Optional numeric length-2 vector for radial limits.
#' @param title Optional title.
#' @param subtitle Optional subtitle.
#' @param legend_title Legend title.
#' @param legend_position Legend position.
#' @param theme Plot theme.
#' @param base_size Base font size.
#'
#' @return A \code{ggplot} object.
#' @export
#'
#' @examples
#' data(pet_maps_dk68)
#' perm_id <- replicate(20, sample.int(nrow(pet_maps_dk68)), simplify = "matrix")
#' res <- pet_corr(
#'   brain = pet_maps_dk68[, c("D1", "DAT", "MOR")],
#'   pet_select = c("D1", "D2", "DAT", "MOR", "5-HT2a"),
#'   perm_id = perm_id,
#'   n_perm = 20
#' )
#' p <- pet_plot_radar(res, metric = "r_obs")
pet_plot_radar <- function(pet_corr_obj,
                           brain_map = NULL,
                           metric = c("r_obs", "minus_log10_p_spin"),
                           pet_order = NULL,
                           line_size = 0.9,
                           point_size = 2.2,
                           fill_alpha = 0.14,
                           palette = NULL,
                           radial_limits = NULL,
                           title = NULL,
                           subtitle = NULL,
                           legend_title = "Brain map",
                           legend_position = "right",
                           theme = c("minimal", "classic", "bw"),
                           base_size = 11) {
  .check_pkg("ggplot2", "PET plotting")
  .pet_plot_result_check(pet_corr_obj)
  metric <- match.arg(metric)
  theme <- match.arg(theme)

  df <- .pet_plot_filter_results(pet_corr_obj, brain_map = brain_map)
  df$value <- if (metric == "r_obs") {
    df$r_obs
  } else {
    -log10(pmax(df$p_spin, .Machine$double.xmin))
  }

  if (!is.null(pet_order)) {
    if (!is.character(pet_order)) {
      .pet_stop("E_PET_PLOT_ORDER", "`pet_order` must be a character vector.")
    }
    missing_pet <- setdiff(unique(df$pet_map), pet_order)
    if (length(missing_pet) > 0) {
      pet_order <- c(pet_order, missing_pet)
    }
  } else {
    meta <- pet_available()
    pet_order <- unique(meta$pet_map[meta$pet_map %in% unique(df$pet_map)])
    if (length(pet_order) == 0) {
      pet_order <- unique(df$pet_map)
    }
  }

  df$pet_map <- factor(df$pet_map, levels = pet_order)
  df <- df[order(df$brain_map, df$pet_map), , drop = FALSE]

  # Close each polygon by appending first point at the end.
  closed <- lapply(split(df, df$brain_map), function(d) {
    d <- d[order(d$pet_map), , drop = FALSE]
    rbind(d, d[1, , drop = FALSE])
  })
  df_closed <- do.call(rbind, closed)
  rownames(df_closed) <- NULL

  p <- ggplot2::ggplot(
    df_closed,
    ggplot2::aes(x = pet_map, y = value, group = brain_map, color = brain_map, fill = brain_map)
  ) +
    ggplot2::geom_polygon(alpha = fill_alpha, linewidth = line_size) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::coord_polar() +
    ggplot2::labs(
      x = NULL,
      y = if (metric == "r_obs") "Observed correlation (r)" else expression(-log[10](p[spin])),
      color = legend_title,
      fill = legend_title,
      title = title %||% "PET correlation radar profile",
      subtitle = subtitle
    ) +
    .pet_plot_theme(theme = theme, base_size = base_size) +
    ggplot2::theme(legend.position = legend_position)

  if (!is.null(radial_limits)) {
    if (!is.numeric(radial_limits) || length(radial_limits) != 2) {
      .pet_stop("E_PET_PLOT_LIMITS", "`radial_limits` must be numeric length 2.")
    }
    p <- p + ggplot2::scale_y_continuous(limits = radial_limits)
  }
  if (!is.null(palette)) {
    p <- p +
      ggplot2::scale_color_manual(values = palette) +
      ggplot2::scale_fill_manual(values = palette)
  }
  p
}
