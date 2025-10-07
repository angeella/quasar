#' @title Closed Testing procedure
#' @description This function applies the closed testing procedure to strong control
#' of the FWER to test the effect of a covariate of interest when multiple
#' quantile regressions are fitted.
#' @usage closedTesting(res, quantiles)
#' @param res data.frame object from \code{\link{waldTest}}, \code{\link{rankTest}}
#' @param quantiles numeric vector that contains the set of quantiles of interest. If \code{NULL}
#' all single hypotheses are considered.
#' @author Angela Andreella
#' @return Data.frame containing the quantile set defined in \code{quantiles} and the related adjusted \eqn{p}-values.
#' @seealso \code{\link{waldTest}}, \code{\link{rankTest}}
#' @export
#' @references
#' Marcus, R., Eric, P., & Gabriel, K. R. (1976). On closed testing procedures with special reference to ordered analysis of variance. Biometrika, 63(3), 655-660.
#'
#' Goeman, J. J., Hemerik, J., & Solari, A. (2021). Only closed testing procedures are admissible for controlling false discovery proportions. The Annals of Statistics, 49(2), 1218-1238.
#'
#' @examples
#' library(quantreg)
#' tau=c(0.1,0.25,0.5,0.75,0.9)
#' x<-rnorm(100)
#' z<-rnorm(100)
#' y<-rnorm(100)
#' mod <- rq(y ~ x + z, tau = tau)
#' res <- waldTest(mod = mod,  X = "x")
#' closedTesting(res, quantiles = c(0.1))


closedTesting <- function(res, quantiles){

  parse_set <- function(s) {
    s <- gsub("[()\\s]", "", as.character(s))
    as.numeric(strsplit(s, ",")[[1]])
  }

  if(!is.null(quantiles)){
  idx <- vapply(res$set, function(s) {
      vals <- parse_set(s)
      all(quantiles %in% vals)
    }, logical(1))


  if(sum(idx) == 0){
    stop("Some quantiles are not the one considered in the model")
  }

  res <- res[idx, ]

  out <- data.frame(p.value.adjusted = max(res$p.value),
                    quantiles.set = paste(quantiles, collapse = ", "))

  }else{

    out <- data.frame(p.value.adjusted = NA,
                      quantiles.set = NA)

    quantiles_single <- res$set[vapply(res$set, function(s) {
      vals <- parse_set(s)
      length(vals) == 1
    }, logical(1))]

    for(i in seq(length(quantiles_single))){

      idx <- vapply(res$set, function(s) {
        vals <- parse_set(s)
        all(quantiles_single[i] %in% vals)
      }, logical(1))

      res_sub <- res[idx, ]

      out[i,] <- c(max(res_sub$p.value),
                   quantiles_single[i])

    }

  }

  return(out)
}
