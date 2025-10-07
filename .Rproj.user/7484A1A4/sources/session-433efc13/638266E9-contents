#' @title Wald-type test for quantile regression
#' @description This function computes the Wald-type test
#' @usage waldTest(mod, X, h = NULL, beta = 0)
#' @param mod \code{rqs} object from \code{rq} function (\code{quantreg} package) that is the
#' quantile models of interest
#' @param X string value that is the covariate of interest
#' @param h bandwidth value, default \code{NULL}
#' @param beta numeric value that is the value of the parameter of interest under the null hypothesis
#' Default 0.
#' @importFrom combinat combn
#' @importFrom stats pchisq
#' @author Angela Andreella
#' @return Data.frame containing the quantile set, statistical test and related \eqn{p}-value
#' @seealso \code{\link{rankTest}}
#' @export

#' @examples
#' library(quantreg)
#' tau=c(0.1,0.25,0.5,0.75,0.9)
#' x<-rnorm(100)
#' z<-rnorm(100)
#' y<-rnorm(100)
#' mod <- rq(y ~ x + z, tau = tau)
#' waldTest(mod = mod,  X = "x")


waldTest <- function(mod, X, h = NULL, beta = 0){

  taus <- mod$tau
  coefest <- mod$coefficients[rownames(mod$coefficients) %in% X,]

  vcov <- estimateCovariance(mod = mod, X = X,test = "wald", h = h)

  design <- mod$x

  tests<-unlist(lapply(1:length(taus),
                       combn,
                       x = taus,
                       simplify = FALSE),
                recursive = FALSE)


  pval<-numeric(0)
  tstat <- numeric(0)

  for (l in 1:length(tests)) {

    coefestsub <- as.matrix(coefest[which(taus %in% unlist(tests[l]))])
    coefestsub <- coefestsub - matrix(beta, nrow = nrow(coefestsub))

    vcovblock <- vcov[c(ncol(design)*which(taus %in% unlist(tests[l]))-1),
                  c(ncol(design)*which(taus %in% unlist(tests[l]))-1)]

    tstat[l] <- t(coefestsub)%*%solve(vcovblock)%*%(coefestsub)

    pval[l] <-1-pchisq(tstat[l],df = length(coefestsub))
  }

  set = gsub(pattern = "c", replacement = "", x = paste0(tests))
  out <- data.frame(p.value = pval,
                    t.stat = tstat,
                    set = set)

  return(out)
}
