#' @title Simulate data
#' @description Simulates covariates \code{(X, Z)} and a response \code{y} drawn from
#' the chosen distribution
#' @usage simulateData(n, beta, gamma, mu, Sigma, sigma.y,
#'   distribution = "normal", df = 5, seed = 1234)
#' @param n Integer. Number of observations.
#' @param beta Numeric scalar. Effect of \code{X}.
#' @param gamma Numeric vector. Effects for \code{Z} (length \code{p - 1}, where \code{p = ncol(Sigma)}).
#' @param mu Numeric scalar. Intercept.
#' @param Sigma Numeric \code{p x p} symmetric positive-definite covariance matrix for \code{(X, Z)}.
#'   The first column corresponds to \code{X}, the remaining columns are \code{Z1, Z2, ...}.
#' @param sigma.y Either a numeric scalar or a one-sided expression/string (e.g., \code{"0.3*abs(X)+0.1"})
#'   defining the scale of \code{y}.
#' @param distribution Character. One of \code{"normal"}, \code{"exponential"}, or \code{"t"}.
#'   This is the distribution of \code{y}.
#' @param df Numeric scalar > 0. Degrees of freedom for t-distribution.
#' @param seed Numeric scalar > 0. Seed for random number generator
#' @author Angela Andreella
#' @return
#' A \code{data.frame} with columns \code{y}, \code{X}, and \code{Z1, ..., Zk}.
#'
#' @examples
#' set.seed(1)
#' p <- 3; Sigma <- diag(p)
#' # Normal
#' dat_n <- simulateData(200, beta=0.5, gamma=c(0.2,-0.1), mu=0,
#'                       Sigma=Sigma, sigma.y=0.5, distribution="normal")
#' # Exponential
#' dat_e <- simulateData(200, beta=0.5, gamma=c(0.2,-0.1), mu=0,
#'                       Sigma=Sigma, sigma.y="0.3*abs(X)+0.1", distribution="exponential")
#' # Student-t
#' dat_t0 <- simulateData(200, beta=0.5, gamma=c(0.2,-0.1), mu=0,
#'                        Sigma=Sigma, sigma.y=0.5, distribution="t", df=7)
#'
#' @importFrom MASS mvrnorm
#' @importFrom stats rnorm rexp rt
#' @export

simulateData <- function(n, beta, gamma, mu, Sigma, sigma.y,
                         distribution = "normal", df = 5, seed = 1234) {

  set.seed(seed)

  p <- ncol(Sigma)

  XZ <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma)

  colnames(XZ) <- c("X", if (p > 1) paste0("Z", seq_len(p - 1)) else character(0))

  eta <- as.numeric(mu + beta * XZ[, 1] + if (p > 1) drop(XZ[, -1, drop = FALSE] %*% gamma) else 0)

  if(is.character(sigma.y)){

    env <- as.data.frame(XZ)
    sigma <- eval(str2lang(sigma.y), envir = env)

  }

  if (is.numeric(sigma.y) && length(sigma.y) == 1L) {
    sigma <- rep(sigma.y, n)
  }

  if (distribution == "exponential") {
    eta_clip <- pmin(pmax(eta, -30), 30)
    rate <- 1 / (exp(eta_clip) * sigma)
    y <- rexp(n, rate = rate) - 1/rate
  }

  if (distribution == "t") {
    y <- eta + sigma * sqrt((df-2)/df) * rt(n = n, df = df)
  }

  if(distribution == "normal"){
    y <- rnorm(n = n, mean = eta, sd = sigma)
  }

  out <- as.data.frame(cbind(y, XZ))
  return(out)
}
