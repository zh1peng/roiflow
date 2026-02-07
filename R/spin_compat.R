# Legacy compatibility wrappers for historical spin-correlation workflow.

.corr_p_matrix <- function(mat, method = c("pearson", "spearman"), use = "pairwise.complete.obs") {
  method <- match.arg(method)
  mat <- .as_numeric_matrix(mat, "df2cor")
  n <- ncol(mat)
  r <- matrix(NA_real_, n, n, dimnames = list(colnames(mat), colnames(mat)))
  p <- matrix(NA_real_, n, n, dimnames = list(colnames(mat), colnames(mat)))

  for (i in seq_len(n)) {
    for (j in i:n) {
      if (i == j) {
        r[i, j] <- 1
        p[i, j] <- NA_real_
      } else {
        xi <- mat[, i]
        xj <- mat[, j]
        keep <- stats::complete.cases(xi, xj)
        if (!any(keep)) {
          r[i, j] <- NA_real_
          p[i, j] <- NA_real_
        } else {
          rr <- stats::cor(xi[keep], xj[keep], method = method, use = use)
          tt <- suppressWarnings(stats::cor.test(xi[keep], xj[keep], method = method))
          r[i, j] <- rr
          p[i, j] <- tt$p.value
        }
        r[j, i] <- r[i, j]
        p[j, i] <- p[i, j]
      }
    }
  }
  list(r = r, P = p)
}

#' Legacy Correlation Matrix Plot with Spin Significance
#'
#' @description
#' Backward-compatible plot helper that mirrors historical behavior.
#' Prefer \code{\link{plot_spin_matrix}} for new code.
#'
#' @param cor.list List containing \code{r} and \code{p.spin}.
#' @param ... Passed to \code{corrplot::corrplot}.
#' @return A corrplot side-effect.
#' @export
cor_plot <- function(cor.list, ...) {
  .check_pkg("corrplot", "legacy spin correlation plots")
  col2use <- grDevices::colorRampPalette(c("steelblue1", "white", "firebrick1"))
  corrplot::corrplot(
    cor.list$r,
    method = "color",
    col = col2use(200),
    addCoef.col = "black",
    tl.col = "black",
    tl.srt = 45,
    p.mat = cor.list$p.spin,
    sig.level = 0.05,
    insig = "blank",
    diag = FALSE,
    number.cex = 0.8,
    cl.ratio = 0.3,
    addgrid.col = "grey",
    ...
  )
}

#' Legacy Correlation Export Table
#'
#' @description
#' Converts legacy \code{cor.list} into an upper-triangle stacked table.
#'
#' @param cor.list List containing \code{r}, \code{P}, and \code{p.spin}.
#' @return Matrix with stacked upper triangles.
#' @export
cor_t2save <- function(cor.list) {
  r2save <- cor.list$r
  p2save <- cor.list$P
  p.spin2save <- cor.list$p.spin
  r2save[lower.tri(r2save)] <- NA
  p2save[lower.tri(p2save)] <- NA
  p.spin2save[lower.tri(p.spin2save)] <- NA
  rbind(r2save, p2save, p.spin2save)
}

#' Legacy Monolithic Spin Workflow Wrapper
#'
#' @description
#' Compatibility wrapper for historical code that combined computation, plotting,
#' and saving in a single function. New code should use:
#' \code{\link{spin_correlation}}, \code{\link{save_spin_results}},
#' \code{\link{plot_spin_matrix}}, and \code{\link{plot_spin_null}}.
#'
#' @param df2cor Numeric matrix/data.frame with maps in columns.
#' @param perm.id Optional permutation matrix.
#' @param nrot Integer permutation count (if \code{perm.id} is \code{NULL}).
#' @param seed Optional seed.
#' @param corr.type Correlation method.
#' @param res_path Optional output folder for legacy file outputs.
#' @param make_plot Logical. Save legacy plot if \code{res_path} is set.
#' @param save_table Logical. Save legacy table if \code{res_path} is set.
#' @param parallel Logical.
#' @param n_cores Integer cores for parallel mode.
#'
#' @return Legacy-style list with \code{r}, \code{P}, \code{p.spin}, and
#' \code{spin_result}.
#' @export
functions_spin_corr <- function(df2cor,
                                perm.id = NULL,
                                nrot = 5000L,
                                seed = 666L,
                                corr.type = c("pearson", "spearman"),
                                res_path = NULL,
                                make_plot = TRUE,
                                save_table = TRUE,
                                parallel = FALSE,
                                n_cores = 2L) {
  warning("`functions_spin_corr()` is deprecated; use `spin_correlation()`.", call. = FALSE)
  corr.type <- match.arg(corr.type)
  df2cor <- .as_numeric_matrix(df2cor, "df2cor")

  if (is.null(perm.id)) {
    perm.id <- spin_generate_permutations(n_perm = nrot, seed = seed)
  }

  spin_res <- spin_correlation(
    x = df2cor,
    y = NULL,
    perm_id = perm.id,
    method = corr.type,
    seed = seed,
    parallel = parallel,
    n_cores = n_cores,
    return_null = FALSE
  )

  cp <- .corr_p_matrix(df2cor, method = corr.type)
  cor.list <- list(
    r = cp$r,
    P = cp$P,
    p.spin = spin_res$p_value,
    spin_result = spin_res
  )

  if (!is.null(res_path) && nzchar(res_path)) {
    dir.create(res_path, recursive = TRUE, showWarnings = FALSE)
    if (isTRUE(make_plot)) {
      grDevices::png(
        file = file.path(res_path, "subtype_vs_disorders_cor_mat_spin_sig.png"),
        width = 7, height = 7, units = "in", res = 300
      )
      cor_plot(cor.list, type = "upper")
      grDevices::dev.off()
    }
    if (isTRUE(save_table)) {
      utils::write.csv(
        cor_t2save(cor.list),
        file.path(res_path, "subtype_vs_disorders_cor_mat_spin.csv"),
        row.names = FALSE
      )
    }
  }

  cor.list
}
