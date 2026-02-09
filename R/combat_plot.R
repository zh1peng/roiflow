# ComBat diagnostic plotting utilities.

.combat_stop <- function(code, msg) {
  stop(sprintf("%s: %s", code, msg), call. = FALSE)
}

.combat_extract_df <- function(x, arg = "data") {
  if (is.data.frame(x)) return(x)
  if (is.list(x) && !is.null(x$data) && is.data.frame(x$data)) return(x$data)
  .combat_stop(
    "E_COMBAT_PLOT_INPUT",
    sprintf("`%s` must be a data.frame or an object containing `$data` as a data.frame.", arg)
  )
}

.combat_infer_study_col <- function(df_before, df_after, study_col) {
  if (!is.null(study_col)) return(study_col)
  if ("Study" %in% names(df_before) && "Study" %in% names(df_after)) return("Study")
  NULL
}

.combat_default_roi <- function(df_before, df_after, site_col, group_col, study_col, roi_col, value_col) {
  long_before <- all(c(roi_col, value_col) %in% names(df_before))
  long_after <- all(c(roi_col, value_col) %in% names(df_after))
  if (long_before && long_after) {
    roi_before <- unique(as.character(df_before[[roi_col]]))
    roi_after <- unique(as.character(df_after[[roi_col]]))
    roi_pool <- intersect(roi_before, roi_after)
    if (length(roi_pool) == 0) {
      .combat_stop(
        "E_COMBAT_PLOT_ROI",
        sprintf("No overlapping ROI names found in `%s` across before/after data.", roi_col)
      )
    }
    preferred <- c("L_bankssts_thickavg", "R_bankssts_thickavg", "banksts")
    hit <- preferred[preferred %in% roi_pool]
    return(if (length(hit) > 0) hit[1] else roi_pool[1])
  }

  reserved <- unique(c(site_col, group_col, study_col))
  cols_before <- setdiff(names(df_before), reserved)
  cols_after <- setdiff(names(df_after), reserved)
  cols <- intersect(cols_before, cols_after)
  if (length(cols) == 0) {
    .combat_stop("E_COMBAT_PLOT_ROI", "No shared ROI-like columns found between before/after data.")
  }

  is_num <- vapply(cols, function(nm) is.numeric(df_before[[nm]]) && is.numeric(df_after[[nm]]), logical(1))
  cols <- cols[is_num]
  if (length(cols) == 0) {
    .combat_stop("E_COMBAT_PLOT_ROI", "No shared numeric ROI columns found between before/after data.")
  }

  preferred <- c("L_bankssts_thickavg", "R_bankssts_thickavg", "banksts")
  hit <- preferred[preferred %in% cols]
  if (length(hit) > 0) return(hit[1])

  roi_like <- cols[grepl("^(L_|R_)|roi|thick|vol|area", cols, ignore.case = TRUE)]
  if (length(roi_like) > 0) return(roi_like[1])
  cols[1]
}

.combat_prepare_stage <- function(df, stage, roi, site_col, group_col, value_col, roi_col, study_col) {
  required <- c(site_col, group_col)
  if (!is.null(study_col)) required <- c(required, study_col)
  .ensure_cols(df, required, context = sprintf("combat_plot_roi(%s)", stage))

  has_long <- all(c(roi_col, value_col) %in% names(df))
  if (has_long) {
    roi_vals <- as.character(df[[roi_col]])
    if (!roi %in% roi_vals) {
      .combat_stop(
        "E_COMBAT_PLOT_ROI",
        sprintf("ROI '%s' was not found in `%s` for `%s` data.", roi, roi_col, stage)
      )
    }
    idx <- roi_vals == roi
    d <- df[idx, , drop = FALSE]
    raw_val <- d[[value_col]]
  } else {
    if (!roi %in% names(df)) {
      .combat_stop(
        "E_COMBAT_PLOT_ROI",
        sprintf(
          "ROI '%s' was not found as a column in `%s` data. If data are long format, set `roi_col` and `value_col`.",
          roi, stage
        )
      )
    }
    d <- df
    raw_val <- d[[roi]]
  }

  if (!is.numeric(raw_val)) {
    val <- suppressWarnings(as.numeric(raw_val))
    bad <- sum(!is.na(raw_val) & is.na(val))
    if (bad > 0) {
      .combat_stop(
        "E_COMBAT_PLOT_VALUE",
        sprintf("ROI values for '%s' in `%s` could not be parsed as numeric.", roi, stage)
      )
    }
  } else {
    val <- raw_val
  }

  out <- data.frame(
    stage = stage,
    site = as.character(d[[site_col]]),
    group = as.character(d[[group_col]]),
    value = as.numeric(val),
    stringsAsFactors = FALSE
  )
  if (!is.null(study_col)) {
    out$study <- as.character(d[[study_col]])
  }

  keep <- stats::complete.cases(out)
  dropped <- sum(!keep)
  if (dropped > 0) {
    warning(
      sprintf("W_COMBAT_PLOT_NA: Dropped %d rows with missing plotting fields (%s stage).", dropped, stage),
      call. = FALSE
    )
  }
  out <- out[keep, , drop = FALSE]
  if (nrow(out) == 0) {
    .combat_stop("E_COMBAT_PLOT_EMPTY", sprintf("No complete plotting rows in `%s` data.", stage))
  }
  out
}

.combat_apply_filters <- function(df_plot, site = NULL, study = NULL, has_study = FALSE) {
  if (!is.null(site)) {
    site <- as.character(site)
    bad_site <- setdiff(site, unique(df_plot$site))
    if (length(bad_site) > 0) {
      .combat_stop(
        "E_COMBAT_PLOT_SITE",
        sprintf("Unknown `site` value(s): %s", paste(bad_site, collapse = ", "))
      )
    }
    df_plot <- df_plot[df_plot$site %in% site, , drop = FALSE]
  }

  if (!is.null(study)) {
    if (!has_study || !"study" %in% names(df_plot)) {
      .combat_stop(
        "E_COMBAT_PLOT_STUDY",
        "`study` filter was provided, but no `study_col` is available in the data."
      )
    }
    study <- as.character(study)
    bad_study <- setdiff(study, unique(df_plot$study))
    if (length(bad_study) > 0) {
      .combat_stop(
        "E_COMBAT_PLOT_STUDY",
        sprintf("Unknown `study` value(s): %s", paste(bad_study, collapse = ", "))
      )
    }
    df_plot <- df_plot[df_plot$study %in% study, , drop = FALSE]
  }

  if (nrow(df_plot) == 0) {
    .combat_stop("E_COMBAT_PLOT_EMPTY", "No rows remain after site/study filtering.")
  }
  df_plot
}

.combat_plot_theme <- function(theme = c("minimal", "classic", "bw"), base_size = 11) {
  theme <- match.arg(theme)
  switch(
    theme,
    minimal = ggplot2::theme_minimal(base_size = base_size),
    classic = ggplot2::theme_classic(base_size = base_size),
    bw = ggplot2::theme_bw(base_size = base_size)
  )
}

.combat_plot_build <- function(df_plot, title, subtitle, xlab, ylab,
                               point_size, point_alpha, box_alpha, box_width,
                               jitter_width, dodge_width,
                               palette, legend_position, theme, base_size) {
  dodge <- ggplot2::position_dodge2(width = dodge_width, preserve = "single")
  jitter_dodge <- ggplot2::position_jitterdodge(
    jitter.width = jitter_width,
    jitter.height = 0,
    dodge.width = dodge_width
  )

  p <- ggplot2::ggplot(
    df_plot,
    ggplot2::aes(x = site_label, y = value, fill = group, color = group)
  ) +
    ggplot2::geom_boxplot(
      position = dodge,
      width = box_width,
      alpha = box_alpha,
      outlier.shape = NA
    ) +
    ggplot2::geom_point(
      position = jitter_dodge,
      size = point_size,
      alpha = point_alpha
    ) +
    ggplot2::labs(
      x = xlab %||% "Site",
      y = ylab,
      title = title,
      subtitle = subtitle,
      fill = "Group",
      color = "Group"
    ) +
    .combat_plot_theme(theme = theme, base_size = base_size) +
    ggplot2::theme(
      legend.position = legend_position,
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)
    )

  if (!is.null(palette)) {
    p <- p +
      ggplot2::scale_fill_manual(values = palette) +
      ggplot2::scale_color_manual(values = palette)
  }
  p
}

.combat_save_plot <- function(plot_obj, split, out_file, width, height, dpi) {
  if (is.null(out_file) || any(!nzchar(out_file))) {
    .combat_stop("E_COMBAT_PLOT_IO", "When `save = TRUE`, provide non-empty `out_file`.")
  }

  if (!isTRUE(split)) {
    ggplot2::ggsave(filename = out_file[1], plot = plot_obj, width = width, height = height, dpi = dpi)
    return(invisible(list(plot = out_file[1])))
  }

  files <- out_file
  if (length(files) == 1) {
    ext <- tools::file_ext(files[1])
    stem <- if (nzchar(ext)) {
      sub(paste0("\\.", ext, "$"), "", files[1])
    } else {
      files[1]
    }
    files <- if (nzchar(ext)) {
      c(paste0(stem, "_before.", ext), paste0(stem, "_after.", ext))
    } else {
      c(paste0(stem, "_before"), paste0(stem, "_after"))
    }
  }
  if (length(files) != 2) {
    .combat_stop("E_COMBAT_PLOT_IO", "For `split = TRUE`, `out_file` must be length 1 or 2.")
  }

  ggplot2::ggsave(filename = files[1], plot = plot_obj$before, width = width, height = height, dpi = dpi)
  ggplot2::ggsave(filename = files[2], plot = plot_obj$after, width = width, height = height, dpi = dpi)
  invisible(list(before = files[1], after = files[2]))
}

#' Plot ROI-Level ComBat Diagnostics (Case vs Control, Before vs After)
#'
#' @description
#' Visual diagnostic for ComBat effects on a single ROI, stratified by
#' site (or study-site) and group (e.g., case vs control), comparing
#' before versus after harmonization.
#'
#' The function accepts wide or long-format tables and can be used directly
#' with data frames produced by \code{\link{combat_fit_df}} /
#' \code{\link{combat_apply_df}} workflows.
#'
#' @param df_before Data frame (or object with \code{$data}) before ComBat.
#' @param df_after Data frame (or object with \code{$data}) after ComBat.
#' @param roi Character scalar ROI name. If \code{NULL}, a package default ROI
#'   is inferred from available columns/labels.
#' @param site_col Character site column name.
#' @param group_col Character group column name.
#' @param value_col Character value column name for long format.
#' @param roi_col Character ROI-name column for long format.
#' @param study_col Optional character study column name. If \code{NULL} and a
#'   \code{"Study"} column exists in both data sets, it is used automatically.
#' @param site Optional character vector to filter sites.
#' @param study Optional character vector to filter studies.
#' @param split Logical. If \code{FALSE} (default), returns one faceted plot
#'   (\code{stage = before/after}). If \code{TRUE}, returns two separate plots
#'   in a \code{combat_roi_plot_split} object.
#' @param ylab Optional y-axis label.
#' @param title Optional plot title.
#' @param subtitle Optional plot subtitle.
#' @param palette Optional named/unnamed vector of group colors.
#' @param point_size Numeric point size.
#' @param point_alpha Numeric point alpha.
#' @param box_alpha Numeric box alpha.
#' @param box_width Numeric box width.
#' @param jitter_width Numeric horizontal jitter width.
#' @param dodge_width Numeric dodge width between groups.
#' @param legend_position Legend position string.
#' @param theme Theme choice: \code{"minimal"}, \code{"classic"}, \code{"bw"}.
#' @param base_size Base text size.
#' @param save Logical. Save output image(s) via \code{ggsave}.
#' @param out_file Output filename. For \code{split = TRUE}, use length-1
#'   (auto-suffixed) or length-2 vector.
#' @param width Plot width for saving (inches).
#' @param height Plot height for saving (inches).
#' @param dpi Plot DPI for saving.
#'
#' @return
#' If \code{split = FALSE}, a \code{ggplot} object.
#'
#' If \code{split = TRUE}, a list with class \code{combat_roi_plot_split}
#' containing \code{$before} and \code{$after} ggplot objects.
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' n <- 80
#' site <- rep(c("S1", "S2"), each = n / 2)
#' group <- rep(rep(c("Control", "Case"), each = n / 4), 2)
#' roi_before <- stats::rnorm(n) + ifelse(site == "S2", 0.8, 0) + ifelse(group == "Case", 0.4, 0)
#' roi_after <- roi_before - ifelse(site == "S2", 0.7, 0) + stats::rnorm(n, sd = 0.05)
#' df_before <- data.frame(Site = site, Group = group, banksts = roi_before)
#' df_after <- data.frame(Site = site, Group = group, banksts = roi_after)
#'
#' p <- combat_plot_roi(df_before, df_after, roi = "banksts")
combat_plot_roi <- function(df_before,
                            df_after,
                            roi = NULL,
                            site_col = "Site",
                            group_col = "Group",
                            value_col = "value",
                            roi_col = "roi",
                            study_col = NULL,
                            site = NULL,
                            study = NULL,
                            split = FALSE,
                            ylab = NULL,
                            title = NULL,
                            subtitle = NULL,
                            palette = NULL,
                            point_size = 1.8,
                            point_alpha = 0.55,
                            box_alpha = 0.35,
                            box_width = 0.7,
                            jitter_width = 0.12,
                            dodge_width = 0.75,
                            legend_position = "right",
                            theme = c("minimal", "classic", "bw"),
                            base_size = 11,
                            save = FALSE,
                            out_file = NULL,
                            width = 10,
                            height = 4.5,
                            dpi = 300) {
  .check_pkg("ggplot2", "ComBat ROI plotting")
  theme <- match.arg(theme)

  df_before <- .combat_extract_df(df_before, arg = "df_before")
  df_after <- .combat_extract_df(df_after, arg = "df_after")
  study_col <- .combat_infer_study_col(df_before, df_after, study_col = study_col)

  if (!is.null(roi)) {
    if (!is.character(roi) || length(roi) != 1 || !nzchar(roi)) {
      .combat_stop(
        "E_COMBAT_PLOT_ROI",
        "`roi` must be a single non-empty ROI name. This helper supports single-ROI diagnostics only."
      )
    }
  } else {
    roi <- .combat_default_roi(
      df_before = df_before,
      df_after = df_after,
      site_col = site_col,
      group_col = group_col,
      study_col = study_col,
      roi_col = roi_col,
      value_col = value_col
    )
    warning(
      sprintf("W_COMBAT_DEFAULT_ROI: `roi` was not provided; using default ROI '%s'.", roi),
      call. = FALSE
    )
  }

  before_plot <- .combat_prepare_stage(
    df = df_before,
    stage = "before",
    roi = roi,
    site_col = site_col,
    group_col = group_col,
    value_col = value_col,
    roi_col = roi_col,
    study_col = study_col
  )
  after_plot <- .combat_prepare_stage(
    df = df_after,
    stage = "after",
    roi = roi,
    site_col = site_col,
    group_col = group_col,
    value_col = value_col,
    roi_col = roi_col,
    study_col = study_col
  )

  df_plot <- rbind(before_plot, after_plot)
  has_study <- !is.null(study_col) && "study" %in% names(df_plot)
  df_plot <- .combat_apply_filters(df_plot, site = site, study = study, has_study = has_study)

  if (has_study) {
    df_plot$site_label <- paste(df_plot$study, df_plot$site, sep = "::")
  } else {
    df_plot$site_label <- df_plot$site
  }

  site_levels <- unique(df_plot$site_label)
  stage_levels <- c("before", "after")
  group_levels <- unique(df_plot$group)
  df_plot$site_label <- factor(df_plot$site_label, levels = site_levels)
  df_plot$stage <- factor(df_plot$stage, levels = stage_levels)
  df_plot$group <- factor(df_plot$group, levels = group_levels)

  if (length(group_levels) != 2) {
    warning(
      sprintf(
        "W_COMBAT_GROUP_LEVELS: `%s` has %d level(s) after filtering; expected 2 for case/control view. Plotting anyway.",
        group_col, length(group_levels)
      ),
      call. = FALSE
    )
  }

  ylab <- ylab %||% roi
  title <- title %||% sprintf("ComBat ROI diagnostic: %s", roi)
  subtitle <- subtitle %||% "Case/control by site, before vs after harmonization"
  xlab <- if (has_study) "Study::Site" else "Site"

  if (!isTRUE(split)) {
    p <- .combat_plot_build(
      df_plot = df_plot,
      title = title,
      subtitle = subtitle,
      xlab = xlab,
      ylab = ylab,
      point_size = point_size,
      point_alpha = point_alpha,
      box_alpha = box_alpha,
      box_width = box_width,
      jitter_width = jitter_width,
      dodge_width = dodge_width,
      palette = palette,
      legend_position = legend_position,
      theme = theme,
      base_size = base_size
    ) +
      ggplot2::facet_wrap(~stage, nrow = 1)

    if (isTRUE(save)) {
      .combat_save_plot(
        plot_obj = p,
        split = FALSE,
        out_file = out_file,
        width = width,
        height = height,
        dpi = dpi
      )
    }
    return(p)
  }

  p_before <- .combat_plot_build(
    df_plot = df_plot[df_plot$stage == "before", , drop = FALSE],
    title = paste0(title, " (before)"),
    subtitle = subtitle,
    xlab = xlab,
    ylab = ylab,
    point_size = point_size,
    point_alpha = point_alpha,
    box_alpha = box_alpha,
    box_width = box_width,
    jitter_width = jitter_width,
    dodge_width = dodge_width,
    palette = palette,
    legend_position = legend_position,
    theme = theme,
    base_size = base_size
  )
  p_after <- .combat_plot_build(
    df_plot = df_plot[df_plot$stage == "after", , drop = FALSE],
    title = paste0(title, " (after)"),
    subtitle = subtitle,
    xlab = xlab,
    ylab = ylab,
    point_size = point_size,
    point_alpha = point_alpha,
    box_alpha = box_alpha,
    box_width = box_width,
    jitter_width = jitter_width,
    dodge_width = dodge_width,
    palette = palette,
    legend_position = legend_position,
    theme = theme,
    base_size = base_size
  )

  out <- structure(list(before = p_before, after = p_after), class = c("combat_roi_plot_split", "list"))

  if (isTRUE(save)) {
    .combat_save_plot(
      plot_obj = out,
      split = TRUE,
      out_file = out_file,
      width = width,
      height = height,
      dpi = dpi
    )
  }

  out
}

#' @export
print.combat_roi_plot_split <- function(x, ...) {
  if (is.list(x) && inherits(x$before, "ggplot") && inherits(x$after, "ggplot")) {
    print(x$before)
    print(x$after)
  }
  invisible(x)
}
