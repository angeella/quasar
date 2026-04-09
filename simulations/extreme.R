# remotes::install_github("angeella/quasar")
# setwd(...)

library(quasar)

# ===============================================
# Simulation Study: Extreme quantiles (Figure 4 )
# ===============================================
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


# Parameters used when B = "distribution" (density-based weighting).
# Keep in one place for transparency / reproducibility.
# NOTE: these parameters are passed to quasar via `error.par`.
parameters_list <- list(mean = 0, sd = 1, df = 5, xi = -1.453, omega = 2, alpha = 2.2)

sim <- expand.grid(
  sim            = seq(nsim),
  distributions  = distributions,
  sigma.y        = sigma.y,
  beta           = beta,
  gamma          = gamma,
  stringsAsFactors = FALSE
)
# We create three copies of the same design grid because we will store results
# for three different B choices in separate data.frames.
sim1 <- sim
sim2 <- sim


# -------------------------
# Pre-allocate output columns
# -------------------------
# We store adjusted p-values for the 3 individual hypotheses (k = 3 quantiles):
# - pvC*: closed testing adjusted p-values
# - pvH*: Holm adjustment of raw p-values
# - pvB*: Bonferroni adjustment of raw p-values
#

# For each run we store p-values for the 3 individual quantile hypotheses
# (one per tau):
#   pvC*  : closed testing adjusted p-values
#   pvRH* : Holm adjustment applied to the raw p-values
#   pvRB* : Bonferroni adjustment applied to the raw p-values
#
# The suffix 1,2,3 corresponds to tau = 0.05, 0.10, 0.15, respectively.
sim$pvC1  <- sim$pvC2  <- sim$pvC3  <- NA_real_
sim$pvRH1 <- sim$pvRH2 <- sim$pvRH3 <- NA_real_
sim$pvRB1 <- sim$pvRB2 <- sim$pvRB3 <- NA_real_

sim1$pvC1  <- sim1$pvC2  <- sim1$pvC3  <- NA_real_
sim1$pvRH1 <- sim1$pvRH2 <- sim1$pvRH3 <- NA_real_
sim1$pvRB1 <- sim1$pvRB2 <- sim1$pvRB3 <- NA_real_

sim2$pvC1  <- sim2$pvC2  <- sim2$pvC3  <- NA_real_
sim2$pvRH1 <- sim2$pvRH2 <- sim2$pvRH3 <- NA_real_
sim2$pvRB1 <- sim2$pvRB2 <- sim2$pvRB3 <- NA_real_

# Main loop

for (i in seq(nrow(sim))) {

  # -----------------------
  # a) Simulate data
  # -----------------------
  dat <- simulateData(
    100,
    beta         = sim$beta[i],
    gamma        = sim$gamma[i],
    mu           = 0.5,
    Sigma        = Sigma,
    sigma.y      = sim$sigma.y[i],
    distribution = sim$distributions[i],
    df           = 5,
    seed         = i,
    xi           = -1.453,
    omega        = 2,
    alpha        = 2.2
  )

  # -----------------------
  # b) Fit quantile regression at extreme taus
  # -----------------------
  mod <- rq(y ~ X + Z1, tau = tau, data = dat)

  # =====================================================================
  # (c1) Closed testing with B = "identity"
  # =====================================================================
  # closedTesting returns:
  #   - out$p.value           : raw marginal p-values for each tau
  #   - out$p.value.adjusted  : closed testing adjusted p-values (strong FWER)
  out <- closedTesting(
    mod,
    X           = "X",
    tau         = tau,
    test        = "rank-score",
    B           = "identity",
    error.distr = sim$distributions[i],
    error.par   = parameters_list
  )

  # Store closed testing adjusted p-values
  sim[i, c("pvC1", "pvC2", "pvC3")] <- out$p.value.adjusted

  # Store Holm adjustment on the 3 raw p-values
  sim[i, c("pvRH1", "pvRH2", "pvRH3")] <- p.adjust(out$p.value, method = "holm")

  # Store Bonferroni adjustment on the 3 raw p-values
  sim[i, c("pvRB1", "pvRB2", "pvRB3")] <- pmin(out$p.value * 3, 1)

  # =====================================================================
  # (c2) Closed testing with B = "inverse diagonal"
  # =====================================================================
  out <- closedTesting(
    mod,
    X           = "X",
    tau         = tau,
    test        = "rank-score",
    B           = "inverse diagonal",
    error.distr = sim$distributions[i],
    error.par   = parameters_list
  )

  sim1[i, c("pvC1", "pvC2", "pvC3")] <- out$p.value.adjusted
  sim1[i, c("pvRH1", "pvRH2", "pvRH3")] <- p.adjust(out$p.value, method = "holm")
  sim1[i, c("pvRB1", "pvRB2", "pvRB3")] <- pmin(out$p.value * 3, 1)

  # =====================================================================
  # (c3) Closed testing with B = "distribution"
  # =====================================================================
  out <- closedTesting(
    mod,
    X           = "X",
    tau         = tau,
    test        = "rank-score",
    B           = "distribution",
    error.distr = sim$distributions[i],
    error.par   = parameters_list
  )

  sim2[i, c("pvC1", "pvC2", "pvC3")] <- out$p.value.adjusted
  sim2[i, c("pvRH1", "pvRH2", "pvRH3")] <- p.adjust(out$p.value, method = "holm")
  sim2[i, c("pvRB1", "pvRB2", "pvRB3")] <- pmin(out$p.value * 3, 1)
}


sim$B <- "identity"
sim1$B <- "inverse diagonal"
sim2$B <- "distribution"

sim <- rbind(sim, sim1, sim2)
