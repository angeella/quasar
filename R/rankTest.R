#' @title Rank-score test for quantile regression
#' @description
#' Performs the rank-score test for the covariate of interest
#' \code{X}, at the quantiles defined in \code{tau}, using a fitted quantile
#' regression model. The test evaluates the null hypothesis that the coefficient of \code{X} is equal to zero
#' against a two-sided alternative, at each specified quantile level.
#' Testing equality to a non-zero value is not yet implemented.
#'
#' @usage rankTest(mod, X, tau = NULL, full = FALSE, h = NULL, alpha = 0.05,
#' eps = c(1e-04,1e-04), A = NULL, error.distr = NULL, error.par = NULL)
#'
#' @param mod An object of class \code{rqs} returned by
#'   \code{\link[quantreg]{rq}}, representing the fitted quantile regression models.
#' @param X A string indicating the covariate of interest.
#' @param tau A numeric vector of quantiles of interest used in \code{mod}.
#' If \code{NULL} (default), all quantiles from the \code{mod} object are considered.
#' @param full Logical. If \code{TRUE}, the function returns the test statistics
#' and corresponding \eqn{p}-values for all intersection hypotheses containing
#' \code{tau}. If \code{FALSE} (default), only the results for the single hypotheses
#' are returned.
#' @param h A numeric value for the bandwidth.
#' @param alpha A numeric value used for bandwidth estimation.
#' Following Koenker (2005), it is typically set equal to the desired significance level.
#' @param eps 2d vector of relative accuracies requested for approximate the distribution of the test statistic
#' by Imhof (1961) "Computing the Distribution of Quadratic Forms in Normal Variables"
#' @param A Optional weighting matrix (or a list of matrices) used in the
#' computation of the test statistic.
#'
#' If \code{NULL} (default), the identity matrix is used for each intersection
#' hypothesis.
#'
#' If a single numeric value is supplied, it is interpreted as a scalar and the
#' weighting matrix is set to \eqn{A = c I_k}.
#'
#' If a list is provided, it must contain one matrix for each tested intersection
#' hypothesis, and each matrix must be square with dimension equal to the size
#' of the corresponding quantile subset.
#'
#' This argument is ignored when \code{error.distr} is specified.
#'
#' @param error.distr A character string specifying the assumed distribution of
#' the regression errors, used to construct the weighting matrix \eqn{A} based on
#' the reciprocal of the error density evaluated at the relevant quantile levels.
#'
#' Allowed values are:
#' \itemize{
#'   \item \code{NULL}: default, no distributional assumption; the argument \code{A} is used instead.
#'   \item \code{"normal"}: normal errors with parameters specified in \code{error.par}.
#'   \item \code{"exponential"}: exponential errors with rate/scale specified in \code{error.par}.
#'   \item \code{"t"}: Student-\eqn{t} errors with degrees of freedom specified in \code{error.par}.
#' }
#'
#' When \code{error.distr} is not \code{NULL}, the argument \code{A} is ignored and the
#' matrix \eqn{A} is automatically set to a diagonal matrix with entries
#' \eqn{1 / f(F^{-1}(\tau_j))}, where \eqn{f} is the density of the specified
#' error distribution and \eqn{F^{-1}} its quantile function.
#'
#' @param error.par A named list of parameters associated with the distribution
#' specified in \code{error.distr}.
#'
#' Required elements depend on the chosen distribution:
#' \itemize{
#'   \item For \code{"normal"}: \code{mean} (default 0), \code{sd} (default 1).
#'   \item For \code{"exponential"}: either \code{rate} or \code{scale} (default \code{rate = 1}).
#'   \item For \code{"t"}: \code{df} (degrees of freedom, required).
#' }
#'
#' @return
#' A \code{data.frame} containing:
#' \itemize{
#'   \item \code{Quantiles.Set}: quantile levels
#'   \item \code{Statistic}: rank-score test statistic
#'   \item \code{p.value}: corresponding unadjusted \eqn{p}-value
#' }
#'
#' @author Angela Andreella
#'
#' @references
#' Koenker, R. (2005). \emph{Quantile Regression}. Cambridge University Press.
#'
#' @details
#' This procedure requires that the covariate of interest \code{X} is either numeric
#' or, if categorical, has at most two levels. Multilevel categorical covariates
#' are not supported and will trigger an error.
#'
#' @importFrom stats pgamma
#' @importFrom utils combn
#' @seealso \code{\link[quantreg]{rq}}, \code{\link[quasar]{waldTest}}
#' @export
#' @examples
#' set.seed(1234)
#' D <- simulateData(n = 100, gamma = 0.5, sigma.y = "1 + 2 * pmax(X, 0)")
#'
#' #Quantile regressions at different levels
#' tau <- c(0.1, 0.25, 0.5, 0.75, 0.9)
#' mod <- quantreg::rq(y ~ X + Z1, tau = tau, data=D)
#'
#' # Rank test
#' rankTest(mod, X = "X")

rankTest <- function(mod, X, tau = NULL, full = FALSE, h = NULL, alpha = 0.05, eps = c(1e-04,1e-04), A = NULL, error.distr = NULL, error.par = NULL){

  assert_intercept_present(mod = mod)
  assert_binary_categorical_X(mod = mod, X = X)

  if(is.null(tau)) tau <- mod$tau
  if (any(!tau %in% mod$tau)) stop("All values in tau must be among the quantiles used in mod.")

  res <- estimateCovariance(mod = mod, X = X, test = "rank", h = h, alpha = alpha)
  S <- res$S
  M <- res$M
  taus <- mod$tau

  tests <- unlist(lapply(1:length(taus),
                       combn,
                       x = taus,
                       simplify = FALSE),
                recursive = FALSE)

  pval<-numeric(0)
  tstat <- numeric(0)

  for (l in 1:length(tests)) {
    this_set <- unlist(tests[l])
    idx      <- which(taus %in% this_set)

    S_sub <- as.matrix(S[idx])
    M_sub <- as.matrix(M[idx, idx])
    Sigma <- M_sub
    k_l   <- length(this_set)

    if (!is.null(error.distr)) {
      Al <- .build_A_from_dist(
        taus   = this_set,
        error.distr = error.distr,
        error.par  = error.par
      )

    } else {

      if (is.null(A)) {
        Al <- diag(1, nrow = k_l, ncol = k_l)

      } else if (is.numeric(A) && length(A) == 1) {

        Al <- diag(A, nrow = k_l, ncol = k_l)

      } else if (is.list(A)) {
        check_dimensions_A(A[[l]], k = k_l)
        Al <- A[[l]]

      } else  if (A == "M_sub"){
        Al <- diag(1/M_sub)
        }else{
        stop("Argument 'A' has unsupported type.", call. = FALSE)
      }
    }

    ASigma   <- Al %*% Sigma
    lambdas  <- svd(ASigma)$d
    lambdas <- lambdas[abs(lambdas)>0.001]
    tstat[l] <- sum(S_sub**2)

    if (length(lambdas) == 1) {
      pval[l] <- 1 - pgamma(tstat[l], shape = 1/2, scale = 2 * lambdas)
    } else {
      pval[l] <- .pImhof(lams = lambdas, x = tstat[l], eps = eps)
    }
  }

  set = gsub(pattern = "c", replacement = "", x = paste0(tests))

  out <- data.frame(Quantiles.Set = set,
                    Statistic = tstat,
                    p.value = pval)

  if(full){
    pat <- paste0("(?<![0-9.])(", paste0(gsub("\\.", "\\\\.", tau), collapse="|"), ")(?![0-9])")
    out <- subset(out, grepl(pat, as.character(set), perl = TRUE))
  }else{
    pat <- paste0("^\\s*(", paste0(gsub("\\.", "\\\\.", tau), collapse="|"), ")\\s*$")
    out <- subset(out, grepl(pat, as.character(set)))
  }

  return(out)
}
