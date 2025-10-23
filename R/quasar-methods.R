#' @title Print and summary methods for quasar objects
#'
#' @description
#' These methods provide basic information about objects of class
#' \code{quasar}, typically returned by the \code{\link{closedTesting}} function.
#'
#' @name quasar-methods
#' @aliases print.quasar summary.quasar
#' @param x An object of class \code{quasar}.
#' @param alpha Significance level.
#' @param ... Additional arguments passed to other methods.
#' @return The input object invisibly.
#' @author Anna Vesely
#'
#' @export
summary.quasar <- function(x, alpha = 0.05, ...) {

  n1 <- sum(x$p.value.adjusted <= alpha, na.rm = TRUE)
  n <- nrow(x)

  cat("Summary of quasar object\n")
  cat("Significant quantiles (level ",
      alpha, "): ", n1, " over ", n, "\n\n", sep = "")

  df <- data.frame(
    Quantile = x$Quantile,
    Coefficient = round(x$Coefficient, 4),
    p.value = round(x$p.value, 4),
    p.value.adjusted = round(x$p.value.adjusted, 4)
  )

  df$Signif <- as.character(cut(
    df$p.value.adjusted,
    breaks = c(-Inf, 0.001, 0.01, 0.05, 0.10, Inf),
    labels = c("***", "**", "*", ".", "")
  ))

  print.data.frame(df, row.names = FALSE)
  invisible(x)
}


#' @rdname quasar-methods
#' @export
print.quasar <- function(x, ...) {
  cat("Object of class quasar\n")
  cat("Number of quantiles:", nrow(x), "\n\n")
  print.data.frame(as.data.frame(x))
  invisible(x)
}


