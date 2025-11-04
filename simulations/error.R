# remotes::install_github("angeella/quasar")
# setwd(...)

library(quasar)
library(quantreg)

# =========================
# Type I error control (Figure 2)
# =========================

# Here, we simulate data under the null hypothesis (beta = 0)
# to check the type I error, i.e., Figure 2 (Section 5)

# Set the values for simulations
nsim  <- 1000    # number of simulations
n     <- 100     # number of observations
tau   <- c(0.1, 0.25, 0.5, 0.75, 0.9)   # quantile levels
Sigma <- matrix(c(1, 0.3, 0.3, 1), 2, 2) # covariance matrix for the covariates
distr <- "normal"                        # distribution of the response
gamma <- 0.5                             # value of the nuisance parameter
beta  <- 0                               # value of the parameter of interest tested

# Create matrices to save the output
res_rank <- matrix(NA, nsim, 31)
res_wald <- matrix(NA, nsim, 31)

# Perform the simulations
for (i in 1:nsim) {

  # Create data
  dat_n <- simulateData(
    n            = n,
    beta         = beta,
    gamma        = gamma,
    mu           = 0.5,
    Sigma        = Sigma,
    sigma.y      = 1,
    distribution = distr,
    seed         = i
  )

  # Fit quantile regression
  mod <- rq(y ~ X + Z1, tau = tau, data = dat_n)

  # Compute rank-score test and save related p-values
  res_rank[i, ] <- rankTest(mod, "X", full = TRUE)$p.value

  # Compute Wald-type test and save related p-values
  res_wald[i, ] <- waldTest(mod, "X", full = TRUE)$p.value
}

# Compute type I errors
TypeIcontrol <- data.frame(
  "Rank_values_1" = apply(res_rank, 2, function(x) mean(x < 0.05)),
  "Wald_values_1" = apply(res_wald, 2, function(x) mean(x < 0.05))
)
