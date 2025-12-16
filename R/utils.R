
#' @importFrom stats as.formula
#' @importFrom stats formula
#' @importFrom stats terms
#' @importFrom pracma  quadgk

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

#' Assert that the fitted model includes an intercept
#'
#' @param mod An object of class "rq" or "rqs".
#' @keywords internal
assert_intercept_present <- function(mod) {
  has_int <- tryCatch({
    if(length(mod$tau)==1){
      rn <- names(stats::coef(mod))
    }else{
      rn <- rownames(stats::coef(mod))
    }
    "(Intercept)" %in% rn
  }, error = function(e) FALSE)

  if (!has_int) {
    stop("This method assumes a model with an intercept. ",
         "Please refit the quantile regression including an intercept (e.g., `y ~ 1 + ...`).",
         call. = FALSE)
  }
}


#' Assert that X is numeric or (if categorical) has at most two levels
#'
#' Works with rq/rqs even when model.frame or mod$x are not stored.
#' @param mod An "rq"/"rqs" object
#' @param X   Character scalar: name of the covariate of interest
#' @keywords internal

assert_binary_categorical_X <- function(mod, X) {
  stopifnot(is.character(X), length(X) == 1L)

  mf <- tryCatch(stats::model.frame(mod),
                 error = function(e) NULL)

  if (!is.null(mf) && X %in% names(mf)) {
    xv <- mf[[X]]
    check_logical <- ((is.factor(xv) && nlevels(xv) > 2L) | (is.character(xv) & length(unique(xv))>2L))
    levels_X <- ifelse(is.factor(xv), nlevels(xv), length(unique(xv)))
    if (check_logical) {
      stop(sprintf(
        "The covariate '%s' is a factor with %d levels. ",
        X, levels_X
      ),
      "This method only supports binary categorical covariates (<= 2 levels). ",
      "Please recode it before fitting the model.",
      call. = FALSE)
    }
    return(invisible(TRUE))
  }


}

#' Assert that the matrix A used in the Imhof approximation has
#' valid dimensions at intersection level
#'
#' @param A The matrix A
#' @param k number of tests
#' @keywords internal
check_dimensions_A <- function(A, k) {

  # UNIVARIATE CASE: A is a vector of length 1
  if (is.vector(A) && length(A) == 1) {
    return(invisible(TRUE))
  }

  # MULTIVARIATE CASE: A must be an k × k matrix
  if (is.matrix(A)) {
    if (nrow(A) == k && ncol(A) == k) {
      return(invisible(TRUE))
    } else {
      stop(sprintf(
        "A is a matrix but not %dx%d (it is %dx%d). Please insert a correct list of matrices A",
        k, k, nrow(A), ncol(A)
      ), call. = FALSE)
    }
  }

  stop("Please insert a correct list of matrices A", call. = FALSE)
}


#' Compute the matrix A from  error densities
#'
#' @param taus The level of the quantiles
#' @param error.distr The distribution of the error terms (normal, exponential and t)
#' @param error.par The related parameters of the specified distribution
#' @keywords internal
#' @importFrom stats qnorm dnorm
#' @importFrom stats qexp dexp qt dt

.build_A_from_dist <- function(taus, error.distr, error.par) {
  error.distr <- match.arg(error.distr, c("normal", "exponential", "t"))

  if (error.distr == "normal") {
    mean <- if (!is.null(error.par$mean)) error.par$mean else 0
    sd   <- if (!is.null(error.par$sd))   error.par$sd   else 1

    q_vals <- qnorm(taus, mean = mean, sd = sd)
    f_vals <- dnorm(q_vals, mean = mean, sd = sd)

  } else if (error.distr == "exponential") {
    if (!is.null(error.par$rate) && !is.null(error.par$scale)) {
      stop("Provide either 'rate' or 'scale', not both, in error.par.", call. = FALSE)
    }

    if (!is.null(error.par$rate)) {
      rate <- error.par$rate
    } else if (!is.null(error.par$scale)) {
      rate <- 1 / error.par$scale
    } else {
      rate <- 1  # default: Exp(rate = 1)
    }

    q_vals <- qexp(taus, rate = rate)
    f_vals <- dexp(q_vals, rate = rate)

  } else if (error.distr == "t") {
    if (is.null(error.par$df)) {
      stop("For 't' errors, you must provide 'df' in error.par.", call. = FALSE)
    }
    df <- error.par$df

    q_vals <- qt(taus, df = df)
    f_vals <- dt(q_vals, df = df)
  }

  if (any(f_vals <= 0)) {
    stop("All resulting densities must be positive.", call. = FALSE)
  }

  A  <- diag(1 / f_vals, nrow = length(taus), ncol = length(taus))
  return(A)
}


#' The Imhof method for x=0 and central variables only.
#' @param lams eigenvalues
#' @param x observed stat test
#' @param eps tolerance vector (one for integrate, one for quadgk)
#' @importFrom stats integrate
#' @importFrom pracma quadgk
#' @importFrom methods is
.pImhof <- function(lams, x, eps = c(1e-10, 1e-10)) {
  lams <- lams[lams != 0]

  if (length(lams) == 0) {
    # degenerate case
    return(list(value = NA_real_, error = NA_real_))
  }

  integrand <- function(u) {          # the Imhof integrand. Domain: 0...Inf
    theta <- 0.5 * colSums(atan(outer(lams,u)))  - 0.5 * as.vector(x) * u
    rho <- exp(colSums(0.25 * log(1 + outer(lams^2,u^2))))
    out <- ifelse(u==0, sum(lams)/2, sin(theta)/(u*rho))
    out
  }
  tr.integrand <- function(v) {       # Transformation of the integrand. Domain: 0...1
    K <- sum(abs(lams))/20            # Scaling constant of the transformation (to make it invariant to the scale of lams)

    0.5 + integrand(-log(1-v)/K) / (pi*K*(1-v))
  }
  rt1 <- max(sqrt(.Machine$double.eps), eps[1])
  rt2 <- max(sqrt(.Machine$double.eps), eps[2])

  res <- try(integrate(tr.integrand, 0, 1, rel.tol = rt1), silent=TRUE)

  if (is(res, "try-error")) {
    res <- try(quadgk(tr.integrand, 0, 1, tol = rt2), silent=TRUE)
  } else {
    out <- res$value
    return(out)
  }
  if (is(res, "try-error")) {
    out <- NA
  } else {
    out <- res
  }
  return(out)
}

