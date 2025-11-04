# remotes::install_github("angeella/quasar")
# setwd("...")

library(quasar)
library(quantreg)

# =========================
# Power analysis (Figure 3)
# =========================
#Here, we simulate data under the alternative hypothesis to #analyze
#the power, i.e., Figure 3 (Section 5)

# ---- Simulation settings ----
nsim  <- 1000                      # number of simulations
n     <- 100                       # sample size
tau   <- c(0.10, 0.25, 0.50, 0.75, 0.90)  # quantile levels
Sigma <- matrix(c(1, 0.3, 0.3, 1), 2, 2)  # covariance matrix for covariates

distributions   <- c("t", "exponential", "normal")  # response distributions
gamma           <- 0.5                               # nuisance parameter
beta            <- c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2)   # parameter of interest
sigma_y_levels  <- c("1+abs(X)", "1")                # sigma.y (kept as character)

# ---- Design grid for simulations ----
sim <- expand.grid(
  sim           = seq_len(nsim),
  distributions = distributions,
  sigma.y       = sigma_y_levels,
  beta          = beta,
  gamma         = gamma,
  stringsAsFactors = FALSE
)

# Scenarios considered in Figure 3 (Section 5):
#  - t / exponential with sigma.y = "1"
#  - normal with sigma.y = "1+abs(X)"
sim <- subset(
  sim,
  (distributions %in% c("t", "exponential") & sigma.y == "1") |
    (distributions == "normal" & sigma.y == "1+abs(X)")
)

# ---- Placeholders for adjusted p-values ----
# Closed testing (rank-score)
sim[c("pvC1","pvC2","pvC3","pvC4","pvC5")] <- NA_real_
# Bonferroni (raw rank-score p-values times number of taus)
sim[c("pvRB1","pvRB2","pvRB3","pvRB4","pvRB5")] <- NA_real_

# ---- Run simulations ----
for (i in seq_len(nrow(sim))) {

  # Simulate data
  dat <- simulateData(
    n            = n,
    beta         = sim$beta[i],
    gamma        = sim$gamma[i],
    mu           = 0.5,
    Sigma        = Sigma,
    sigma.y      = sim$sigma.y[i],           # already character
    distribution = sim$distributions[i],
    df           = 5,
    seed         = i
  )

  # Fit quantile regressions
  mod <- rq(y ~ X + Z1, tau = tau, data = dat)

  # Closed testing (rank-score) to control FWER
  out <- closedTesting(mod, X = "X", test = "rank-score")

  # Store adjusted p-values (closed testing)
  sim[i, c("pvC1","pvC2","pvC3","pvC4","pvC5")] <- as.numeric(out$p.value.adjusted)

  # Store Bonferroni-adjusted p-values (cap at 1)
  bonf <- pmin(as.numeric(out$p.value) * length(tau), 1)
  sim[i, c("pvRB1","pvRB2","pvRB3","pvRB4","pvRB5")] <- bonf
}

# ---- Summaries: proportion of p <= 0.05 by scenario ----
# Select the p-value columns (closed testing + Bonferroni)
pv_cols <- grep("^pv(C|RB)[1-5]$", names(sim), value = TRUE)

# Convert to 0/1 indicators for significance and aggregate by groups
X <- data.frame(lapply(sim[pv_cols], function(z) as.numeric(z <= 0.05)))

res <- aggregate(
  X,
  by = list(
    distributions = sim$distributions,
    sigma.y       = sim$sigma.y,
    beta          = sim$beta
  ),
  FUN = mean,
  na.rm = TRUE
)

# Optional: order rows by grouping variables
res <- res[order(res$distributions, res$sigma.y, res$beta), ]

# 'res' contains the empirical power (proportion of p <= 0.05) for each method and scenario.
