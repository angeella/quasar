#' @importFrom quantreg rq
#' @importFrom quantreg bandwidth.rq
#' @importFrom quantreg rq.fit


estimateDensity <- function(mod, tau, X, test = "rank", h = NULL){

  if(test == "rank"){
     formula <- make_h0_formula(mod, X = X)
     mod<-rq(formula, tau=tau, data = mod$model)
  }


  design <- mod$x
  y <- mod$y
  p <- ncol(design)
  n <- nrow(design)
  eps <- .Machine$double.eps^(1/2)
  resid <- mod$residuals

  pz <- sum(abs(resid) < eps)
  if(is.null(h) | identical(h, "Hall-Sheather")){
    h <- bandwidth.rq(tau, n, hs = TRUE)
  }
  if(identical(h, "Bofinger")){
    h <- bandwidth.rq(tau, n, hs = FALSE)
  }

  while ((tau - h < 0) || (tau + h > 1)) h <- h/2

  bhi <- rq.fit(design, y, tau = tau + h)$coef
  blo <- rq.fit(design, y, tau = tau - h)$coef
  dyhat <- design %*% (bhi - blo)

  f_quant_tau <- pmax(0.1, (2 * h)/(dyhat - eps))


    return(f_quant_tau)

}
