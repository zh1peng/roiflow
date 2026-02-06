# Summary helpers

.summarize_groups <- function(df, vars, group_col) {
  stopifnot(is.data.frame(df))
  if (!group_col %in% names(df)) stop("group_col not found in df.", call. = FALSE)
  vars <- intersect(vars, names(df))
  if (length(vars) == 0) return(data.frame())

  g <- factor(as.character(df[[group_col]]))
  g_levels <- levels(g)

  out <- vector("list", length(vars) * length(g_levels))
  k <- 1
  for (v in vars) {
    x <- df[[v]]
    if (!is.numeric(x)) {
      suppressWarnings(x <- as.numeric(x))
    }
    for (lv in g_levels) {
      idx <- which(g == lv)
      xv <- x[idx]
      n <- sum(!is.na(xv))
      m <- mean(xv, na.rm = TRUE)
      s <- stats::sd(xv, na.rm = TRUE)
      se <- if (n > 0) s / sqrt(n) else NA_real_
      out[[k]] <- data.frame(
        var = v,
        mean = m,
        sd = s,
        se = se,
        n = n,
        stringsAsFactors = FALSE
      )
      out[[k]][[group_col]] <- lv
      k <- k + 1
    }
  }

  res <- do.call(rbind, out)
  # stable column order: group_col first
  res <- res[, c(group_col, "var", "mean", "sd", "se", "n"), drop = FALSE]
  rownames(res) <- NULL
  res
}
