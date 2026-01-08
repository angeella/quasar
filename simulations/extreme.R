# remotes::install_github("angeella/quasar")
# setwd(...)

library(quasar)

# =========================
# Simulation Study: Extreme quantiles (Figure 4 )
# =========================
# This script reproduces the simulation study focusing on extreme quantile levels
# tau = {0.05, 0.10, 0.15}, comparing closed testing versus Holm and Bonferroni,
# and evaluating different weighting matrices B under several DGPs.


# Set the values for simulations
nsim  <- 1000    # number of simulations
n     <- 100     # number of observations
Sigma <- matrix(c(1, 0.3, 0.3, 1), 2, 2) # covariance matrix for the covariates
tau <- c(0.05, 0.10, 0.15) # Quantile levels (extreme quantiles)

# Data-generating mechanisms
distributions <- c("normal", "t", "skew-normal")
sigma_y_specs <- c("1", "1+abs(X)", "1+abs(Z1)")
beta_grid     <- c(0.2, 0.4, 0.6, 0.8, 1.2)
gamma_grid    <- c(0.5)

# Weighting matrices for the generalized rank-score test inside closed testing
B_specs <- c("identity", "inverse diagonal", "distribution")

# Parameters used when B = "distribution" (density-based weighting).
# Keep in one place for transparency / reproducibility.
# NOTE: these parameters are passed to quasar via `error.par`.
error_par <- list(mean = 0, sd = 1, df = 5, xi = -1.453, omega = 2, alpha = 2.2)

# Simulation design
sim <- expand.grid(
  sim           = seq_len(nsim),
  distribution  = distributions,
  sigma_y       = sigma_y_specs,
  beta          = beta_grid,
  gamma         = gamma_grid,
  B             = B_specs,
  stringsAsFactors = FALSE
)

# -------------------------
# Pre-allocate output columns
# -------------------------
# We store adjusted p-values for the 3 individual hypotheses (k = 3 quantiles):
# - pvC*: closed testing adjusted p-values
# - pvH*: Holm adjustment of raw p-values
# - pvB*: Bonferroni adjustment of raw p-values
#

init_cols <- function(prefix) {
  sim[[paste0(prefix, "1")]] <- NA_real_
  sim[[paste0(prefix, "2")]] <- NA_real_
  sim[[paste0(prefix, "3")]] <- NA_real_
}

init_cols("pvC")
init_cols("pvH")
init_cols("pvB")


# Main loop

for (i in seq_len(nrow(sim))) {

  # -----------------------
  # 1) Simulate data
  # -----------------------
  dat <- simulateData(
    n            = n,
    beta         = sim$beta[i],
    gamma        = sim$gamma[i],
    mu           = 0.5,
    Sigma        = Sigma,
    sigma.y      = sim$sigma_y[i],
    distribution = sim$distribution[i],
    # distribution-specific parameters (safe to pass even if unused)
    df           = 5,
    xi           = -1.453,
    omega        = 2,
    alpha        = 2.2,
    seed         = i
  )

  # -----------------------
  # 2) Fit quantile regression at extreme taus
  # -----------------------
  mod <- rq(y ~ X + Z1, tau = tau, data = dat)

  # -----------------------
  # 3) Closed testing: generalized rank-score test with chosen B
  # -----------------------
  # When B = "distribution", quasar uses error.distr + error.par to build weights.
  out <- closedTesting(
    mod,
    X          = "X",
    tau        = tau,
    test       = "rank-score",
    B          = sim$B[i],
    error.distr= sim$distribution[i],
    error.par  = error_par
  )

  # -----------------------
  # 4) Store results
  # -----------------------
  # Closed testing adjusted p-values
  sim[i, c("pvC1", "pvC2", "pvC3")] <- out$p.value.adjusted

  # Holm and Bonferroni based on raw p-values for the 3 individual hypotheses
  sim[i, c("pvH1", "pvH2", "pvH3")] <- p.adjust(out$p.value, method = "holm")
  sim[i, c("pvB1", "pvB2", "pvB3")] <- pmin(out$p.value * length(tau), 1)

}


