#' @title simulate data
#' @description This function simulate data
#' @usage simulateData(n, beta, gamma, mu, Sigma, sigma.y)
#' @param n number of observations
#' @param beta value of the parameter of interest
#' @param gamma vector of values of the nuisance parameters
#' @param mu vector of means for covariates
#' @param Sigma covariance matrix for covariates
#' @param sigma.y vector of standard deviations for the response variable
#' @importFrom MASS mvrnorm
#' @importFrom stats rnorm
#' @author Angela Andreella
#' @return matrix
#' @export




simulateData <- function(n, beta, gamma, mu, Sigma, sigma.y){

  X<-mvrnorm(n,
             mu,
             Sigma)


  y <- rnorm(n, beta*X[,1] + X[,-1] %*% gamma, sigma.y)

  sim <- data.frame(X = X,
                    y = y)

  return(sim)
}
