
#' @importFrom stats as.formula
#' @importFrom stats formula
#' @importFrom stats terms

make_h0_formula <- function(mod, X) {
  f  <- formula(mod)
  tt <- terms(f)
  tl <- attr(tt, "term.labels")
  fac <- attr(tt, "factors")
  intc <- attr(tt, "intercept")

  if (length(tl) == 0L) {
    rhs <- if (intc == 1) "1" else "0"
    resp <- as.character(f)[2L]
    return(as.formula(paste(resp, "~", rhs)))
  }

  term_vars <- lapply(seq_along(tl), function(j) {
    rownames(fac)[fac[, j] != 0]
  })

  drop_idx <- which(vapply(term_vars, function(v) any(v %in% X), logical(1)))
  X_idx <- setdiff(seq_along(tl), drop_idx)

  rhs_terms <- tl[X_idx]
  rhs <- paste(c(if (intc == 1) "1" else "0", rhs_terms), collapse = " + ")
  rhs <- sub("^1 \\+ \\s*$", "1", rhs)
  rhs <- sub("^0 \\+ \\s*$", "0", rhs)

  resp <- as.character(f)[2L]
  as.formula(paste(resp, "~", if (nzchar(rhs)) rhs else if (intc == 1) "1" else "0"))
}
