# Spin-correlation plotting helpers.

#' Plot Spin Correlation Matrix
#'
#' @description
#' Heatmap plot for observed spin correlations or spin p-values.
#'
#' @param result A \code{spin_correlation_result} object.
#' @param value One of \code{"observed"} or \code{"p_value"}.
#' @param triangle \code{"full"}, \code{"upper"}, or \code{"lower"}.
#' @param title Optional title.
#' @param midpoint Midpoint for diverging scales (\code{observed} only).
#' @param limits Optional numeric limits (\code{observed} only).
#' @param plot_spec Optional list with \code{theme} (\code{"minimal"},
#' \code{"classic"}, \code{"bw"}), \code{base_size}, and \code{legend_position}.
#'
#' @return A \code{ggplot} object.
#' @export
plot_spin_matrix <- function(result,
                             value = c("observed", "p_value"),
                             triangle = c("full", "upper", "lower"),
                             title = NULL,
                             midpoint = 0,
                             limits = NULL,
                             plot_spec = list()) {
  .check_pkg("ggplot2", "spin plotting")
  value <- match.arg(value)
  triangle <- match.arg(triangle)

  if (!inherits(result, "spin_correlation_result")) {
    .spin_stop("E_SPIN_RESULT", "`result` must be a `spin_correlation_result`.")
  }
  mat <- if (value == "observed") result$observed else result$p_value
  if (!is.matrix(mat)) {
    .spin_stop("E_SPIN_RESULT", sprintf("`%s` matrix not found in result.", value))
  }

  rn <- rownames(mat) %||% paste0("x_", seq_len(nrow(mat)))
  cn <- colnames(mat) %||% paste0("y_", seq_len(ncol(mat)))
  df <- expand.grid(i = seq_len(nrow(mat)), j = seq_len(ncol(mat)), KEEP.OUT.ATTRS = FALSE)
  df$x_map <- factor(rn[df$i], levels = rev(rn))
  df$y_map <- factor(cn[df$j], levels = cn)
  df$value <- as.numeric(mat[cbind(df$i, df$j)])

  if (triangle != "full" && nrow(mat) == ncol(mat) && identical(rn, cn)) {
    if (triangle == "upper") {
      df <- df[df$i <= df$j, , drop = FALSE]
    } else {
      df <- df[df$i >= df$j, , drop = FALSE]
    }
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = y_map, y = x_map, fill = value)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.25) +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      fill = if (value == "observed") "r" else "p_spin",
      title = title %||% if (value == "observed") "Spin correlation matrix" else "Spin p-value matrix"
    )

  if (value == "observed") {
    if (is.null(limits)) {
      mx <- max(abs(df$value), na.rm = TRUE)
      if (!is.finite(mx) || mx == 0) mx <- 1
      limits <- c(-mx, mx)
    }
    p <- p + ggplot2::scale_fill_gradient2(
      low = "steelblue1", mid = "white", high = "firebrick1",
      midpoint = midpoint, limits = limits, na.value = "grey90"
    )
  } else {
    p <- p + ggplot2::scale_fill_gradient(low = "firebrick1", high = "white", limits = c(0, 1), na.value = "grey90")
  }

  theme_name <- plot_spec$theme %||% "minimal"
  base_size <- plot_spec$base_size %||% 11
  p <- p + switch(
    theme_name,
    minimal = ggplot2::theme_minimal(base_size = base_size),
    classic = ggplot2::theme_classic(base_size = base_size),
    bw = ggplot2::theme_bw(base_size = base_size),
    ggplot2::theme_minimal(base_size = base_size)
  )
  if (!is.null(plot_spec$legend_position)) {
    p <- p + ggplot2::theme(legend.position = plot_spec$legend_position)
  }
  p
}

#' Plot Null Distribution for a Selected Map Pair
#'
#' @description
#' Histogram + density plot of permutation null distributions for one map pair.
#' Useful for two-map runs or to inspect one pair from multi-map runs.
#'
#' @param result A \code{spin_correlation_result} object with nulls.
#' @param x_map Optional x-map name. Default: first map.
#' @param y_map Optional y-map name. Default: first map.
#' @param bins Histogram bins.
#' @param title Optional title.
#' @param plot_spec Optional list: theme/base_size/legend_position.
#'
#' @return A \code{ggplot} object.
#' @export
plot_spin_null <- function(result,
                           x_map = NULL,
                           y_map = NULL,
                           bins = 40,
                           title = NULL,
                           plot_spec = list()) {
  .check_pkg("ggplot2", "spin plotting")
  if (!inherits(result, "spin_correlation_result")) {
    .spin_stop("E_SPIN_RESULT", "`result` must be a `spin_correlation_result`.")
  }
  if (is.null(result$null_xy) || is.null(result$null_yx)) {
    .spin_stop("E_SPIN_NULL", "Null distributions are missing (`return_null = FALSE`).")
  }

  x_names <- rownames(result$observed) %||% paste0("x_", seq_len(nrow(result$observed)))
  y_names <- colnames(result$observed) %||% paste0("y_", seq_len(ncol(result$observed)))
  if (is.null(x_map)) x_map <- x_names[1]
  if (is.null(y_map)) y_map <- y_names[1]

  i <- match(x_map, x_names)
  j <- match(y_map, y_names)
  if (is.na(i) || is.na(j)) {
    .spin_stop("E_SPIN_PAIR", "Requested x_map/y_map not found.")
  }

  d <- data.frame(
    value = c(result$null_xy[, i, j], result$null_yx[, i, j]),
    source = rep(c("permute_x", "permute_y"), each = dim(result$null_xy)[1]),
    stringsAsFactors = FALSE
  )
  obs <- result$observed[i, j]

  p <- ggplot2::ggplot(d, ggplot2::aes(x = value, fill = source, color = source)) +
    ggplot2::geom_histogram(position = "identity", alpha = 0.25, bins = bins) +
    ggplot2::geom_density(alpha = 0.15, linewidth = 0.8) +
    ggplot2::geom_vline(xintercept = obs, linewidth = 0.8, linetype = "longdash") +
    ggplot2::labs(
      x = "Null correlation",
      y = "Count / density",
      title = title %||% sprintf("Spin null: %s vs %s", x_map, y_map),
      subtitle = sprintf("Observed r = %.3f; p_spin = %.4f", obs, result$p_value[i, j]),
      fill = "Null source",
      color = "Null source"
    )

  theme_name <- plot_spec$theme %||% "minimal"
  base_size <- plot_spec$base_size %||% 11
  p <- p + switch(
    theme_name,
    minimal = ggplot2::theme_minimal(base_size = base_size),
    classic = ggplot2::theme_classic(base_size = base_size),
    bw = ggplot2::theme_bw(base_size = base_size),
    ggplot2::theme_minimal(base_size = base_size)
  )
  if (!is.null(plot_spec$legend_position)) {
    p <- p + ggplot2::theme(legend.position = plot_spec$legend_position)
  }
  p
}
