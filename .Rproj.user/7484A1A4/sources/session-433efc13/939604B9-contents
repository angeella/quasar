#' @title Rank score test for quantile regression
#' @description This function computes the rank score test
#' @usage rankTest(mod, X, h = NULL)
#' @param mod \code{rqs} object from \code{rq} function (\code{quantreg} package) that is the
#' quantile models of interest
#' @param X string value that is the covariate of interest
#' @param h bandwidth value, default \code{NULL}
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
#' rankTest(mod = mod,  X = "x")

rankTest <- function(mod, X, h = NULL){

  res <- estimateCovariance(mod = mod, X = X, test = "rank", h = h)
  S <- res$S
  M <- res$M
  taus <- mod$tau

  tests<-unlist(lapply(1:length(taus),
                       combn,
                       x = taus,
                       simplify = FALSE),
                recursive = FALSE)

  pval<-numeric(0)
  tstat <- numeric(0)

  for (l in 1:length(tests)) {
    S_sub <- as.matrix(S[which(taus %in% unlist(tests[l]))])

    M_sub <- as.matrix(M[which(taus %in% unlist(tests[l])),
                         which(taus %in% unlist(tests[l]))])

    tstat[l] <- (t(S_sub)%*%solve(M_sub)%*%S_sub)

    pval[l] <- 1-pchisq(tstat[l],df = length(S_sub))
  }

  set = gsub(pattern = "c", replacement = "", x = paste0(tests))
  out <- data.frame(p.value = pval,
             t.stat = tstat,
             set = set)

  return(out)
}
