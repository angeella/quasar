#' @title Closed testing procedure
#'
#' @description
#' Applies the closed testing procedure to strongly control the familywise error rate (FWER)
#' when testing the effect of a covariate of interest across multiple quantile regression models.
#'
#' @usage closedTesting(mod, X, tau = NULL, test = "rank-score", ...)
#'
#' @param mod An object of class \code{"rqs"} returned by
#'   \code{\link[quantreg]{rq}}, representing the fitted quantile regression models.
#' @param tau A numeric vector of quantiles of interest used in \code{mod}.
#'   If \code{NULL} (default), all quantiles from the \code{mod} object are considered.
#' @param X A string indicating the covariate of interest.
#' @param test Character. Type of test to be used. Options are
#'   \code{"wald"} and \code{"rank-score"}. Default is \code{"rank-score"}.
#' @param ... additional arguments, see \code{\link[quasar]{waldTest}}, \code{\link[quasar]{rankTest}}
#' @return
#' A \code{data.frame} containing:
#' \itemize{
#'   \item \code{Quantile}: quantile considered
#'   \item \code{Coefficient}: estimated coefficient
#'   \item \code{Statistic}: test statistic
#'   \item \code{p.value}: raw \eqn{p}-value
#'   \item \code{p.value.adjusted}: adjusted \eqn{p}-value from the closed testing procedure
#' }
#'
#' @seealso
#' \code{\link[quasar]{waldTest}}, \code{\link[quasar]{rankTest}}, \code{\link[quantreg]{rq}}
#'
#' @references
#' Marcus, R., Eric, P., & Gabriel, K. R. (1976).
#' On closed testing procedures with special reference to ordered analysis of variance.
#' \emph{Biometrika}, 63(3), 655--660.
#'
#' Goeman, J. J., Hemerik, J., & Solari, A. (2021).
#' Only closed testing procedures are admissible for controlling false discovery proportions.
#' \emph{The Annals of Statistics}, 49(2), 1218--1238.
#'
#' @author Angela Andreella
#'
#' @details
#' This procedure requires that the covariate of interest \code{X} is either numeric
#' or, if categorical, has at most two levels. Multilevel categorical covariates
#' are not supported and will trigger an error.
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' library(quantreg)
#' tau <- c(0.1, 0.25, 0.5, 0.75, 0.9)
#' n <- 100
#' x <- rnorm(n)
#' z <- rnorm(n)
#' # signal on x to make the example informative
#' y <- 0.5 * x + 0.2 * z + rnorm(n)
#' mod <- rq(y ~ x + z, tau = tau)
#' closedTesting(mod, X = "x")                     # default rank-score
#' closedTesting(mod, X = "x", test = "wald")      # Wald version

closedTesting <- function(mod, X, tau = NULL, test = "rank-score", ...){

  if(is.null(tau)){
    tau <- mod$tau
  }

  if(sum(!(tau %in% mod$tau))!=0){
    stop("The quantiles specified are not the ones used in the quantile regression. Please specify a vector of proper quantiles.")
  }


  parse_set <- function(s) {
    s <- gsub("[()\\s]", "", as.character(s))
    as.numeric(strsplit(s, ",")[[1]])
  }

  if(test == "rank-score"){
    res <- rankTest(mod = mod, X = X, tau = mod$tau, full = TRUE, ...)
  }
  if(test == "wald"){
    res <- waldTest(mod = mod, X = X, tau = mod$tau, full = TRUE, ...)
  }

    out <- data.frame(p.value.adjusted = NA,
                      quantiles.set = NA)

    quantiles_single <- as.character(mod$tau)

    for(i in seq(length(quantiles_single))){

      idx <- vapply(res$Quantiles.Set, function(s) {
        vals <- parse_set(s)
        all(quantiles_single[i] %in% vals)
      }, logical(1))

      res_sub <- res[idx, ]

      out[i,] <- c(max(res_sub$p.value),
                   quantiles_single[i])


  }

  if(length(mod$tau)==1){
    out <- data.frame(Quantile = out$quantiles.set,
                      Coefficients = mod$coefficients[names(mod$model) == X],
                      Statistic = res$Statistic[1:length(mod$tau)],
                      p.value = round(res$p.value[1:length(mod$tau)],7),
                      p.value.adjusted = round(as.numeric(out$p.value.adjusted),7))

  }else{
    out <- data.frame(Quantile = out$quantiles.set,
                      Coefficients = mod$coefficients[names(mod$model) == X,],
                      Statistic = res$Statistic[1:length(mod$tau)],
                      p.value = round(res$p.value[1:length(mod$tau)],7),
                      p.value.adjusted = round(as.numeric(out$p.value.adjusted),7))

    rownames(out) <- NULL


    out <- subset(out, out$Quantile %in% tau)
  }


  return(out)
}
