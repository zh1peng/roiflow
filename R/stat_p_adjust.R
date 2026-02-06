# Multiple-testing utilities

.apply_p_adjust <- function(tbl, p_col = "p_value", spec = analysis_spec()) {
  if (!p_col %in% names(tbl)) {
    return(tbl)
  }
  p_adj_spec <- spec$p_adjust %||% list(method = "fdr", scope = "within_call")
  method <- p_adj_spec$method %||% "fdr"
  scope <- p_adj_spec$scope %||% "within_call"
  scope_id <- p_adj_spec$family_id %||% NA_character_
  scope_desc <- p_adj_spec$family_desc %||% NA_character_

  if (is.null(method) || identical(scope, "none")) {
    tbl$p_adj <- tbl[[p_col]]
  } else if (is.character(scope) && scope %in% names(tbl)) {
    tbl$p_adj <- ave(tbl[[p_col]], tbl[[scope]], FUN = function(x) stats::p.adjust(x, method))
    scope_id <- scope
    scope_desc <- paste0("within_", scope)
  } else {
    tbl$p_adj <- stats::p.adjust(tbl[[p_col]], method)
    if (is.na(scope_id)) scope_id <- "within_call"
    if (is.na(scope_desc)) scope_desc <- "within_call"
  }

  tbl$p_adj_method <- method
  tbl$p_adj_scope_id <- scope_id
  tbl$p_adj_scope_desc <- scope_desc
  tbl
}

