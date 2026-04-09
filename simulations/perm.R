
# =========================
# Simulation Study: Comparison with permutation-based approach by (Mrkvicka et. al 2026 (Figure 4)
# =========================
# This script reproduces the simulation study over 0.1, 0.25, 0.5, 0.75, 0.9 quantiles
# comparing closed testing versus the permutation-based approach by Mrkvicka et. al 2026,


library(quasar)
library(GET)
library(r41sqrt10)
library(tidyverse)

# Set the values for simulations
nsim <- 1000
rho <- 0.3
p <- 2
Sigma <-  (1 - rho) * diag(p) + rho * matrix(1, p, p)

sigma.y <- c("1", "1+abs(X)", "1+abs(Z1)")
beta <- c(0, 0.4, 0.6)
tau <- c(0.1, 0.25, 0.5, 0.75, 0.9)

# Simulation design
sim <- expand.grid(sim = seq(nsim),
                   sigma.y = sigma.y,
                   beta = beta,
                   stringsAsFactors = F)


for(i in seq(nrow(sim))){
  # -----------------------
  # 1) Simulate data
  # -----------------------
  dat <- simulateData(100,
                      beta=sim$beta[i],
                      gamma=0.5,
                      mu=0.5,
                      Sigma=Sigma,
                      sigma.y=sim$sigma.y[i],
                      distribution="normal",
                      seed = i)
  print(i)

  # -----------------------
  # 2) Permutation-based approach
  # -----------------------
  mod_perm <- global_rq(1000, formula.full = y ~ X + Z1,
                        formula.reduced = y ~ Z1,
                        tau = tau, data=dat,
                        permutationstrategy = "remove location scale",
                        GET.args = list(typeone = "fwer"))

  # Logical values (Accept/Rejection)
  sim[i,c("pvPerm1", "pvPerm2", "pvPerm3", "pvPerm4", "pvPerm5")]<-mod_perm$obs<mod_perm$lo | mod_perm$obs >mod_perm$hi

  # -----------------------
  # 3) Fit quantile regression
  # -----------------------
  mod <- rq(y ~ X + Z1, tau = tau, data = dat)

  # -----------------------
  # 4) Closed testing: generalized rank-score test
  # -----------------------
  out <- closedTesting(mod, X = "X", test = "rank-score", tau = tau, B = "identity")

  # Closed testing adjusted p-values
  sim[i, c("pvC1", "pvC2", "pvC3", "pvC4", "pvC5")] <- out$p.value.adjusted


}

#Compute FWER

sim %>%
  filter(beta == 0) %>%
  mutate(
    fwerC = apply(cbind(pvC1, pvC2, pvC3, pvC4, pvC5) <= 0.05, 1, any),
    fwerPerm = apply(cbind(pvPerm1, pvPerm2, pvPerm3, pvPerm4, pvPerm5), 1, any)
  ) %>%
  group_by(sigma.y) %>%
  summarise(
    fwerC = mean(fwerC),
    fwerPerm = mean(fwerPerm),
    .groups = "drop"
  )

#Figure 2

q_map <- c("1" = 0.10,
           "2" = 0.25,
           "3" = 0.50,
           "4" = 0.75,
           "5" = 0.90)

db <- sim %>%
  filter(beta != 0) %>%
  group_by(sigma.y, beta) %>%
  summarize(pvC1 = mean(pvC1 <= 0.05),
            pvC2 = mean(pvC2 <= 0.05),
            pvC3 = mean(pvC3 <= 0.05),
            pvC4 = mean(pvC4 <= 0.05),
            pvC5 = mean(pvC5 <= 0.05),
            pvPerm1 = mean(pvPerm1),
            pvPerm2 = mean(pvPerm2),
            pvPerm3 = mean(pvPerm3),
            pvPerm4 = mean(pvPerm4),
            pvPerm5 = mean(pvPerm5),
            .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(
      starts_with("pvC"),
      starts_with("pvPerm")
    ),
    names_to = c("Methods", "q_id"),
    names_pattern = "^pv([A-Za-z]+)([1-5])$",
    values_to = "value"
  ) %>%
  mutate(
    Methods = case_when(
      Methods == "C"   ~ "Closed Testing",
      Methods == "Perm"   ~ "Permutation-based approach (Mrkvicka et. al 2026)"
    ),
    quantile = q_map[q_id],
    Beta     = as.numeric(beta)
  )

beta_vals <- sort(unique(db$beta))
beta_map  <- setNames(paste0("beta == ", beta_vals), as.character(beta_vals))

lab_sigma <- c(
  "1" = "sigma^2 == 1",
  "1+abs(X)" = "sigma^2 == 1 + abs(X)",
  "1+abs(Z1)" = "sigma^2 == 1 + abs(Z)"
)

db %>%
  mutate(
    quantile = as.character(quantile),
    beta = as.character(beta)
  ) %>%
  ggplot(aes(x = quantile, y = value)) +
  geom_point(aes(x = quantile, y = value, col = Methods), size = 2) +
  geom_line(aes(x = quantile, y = value, group = Methods, col = Methods), size = 1)+
  facet_grid(
    beta ~ sigma.y,
    labeller = labeller(
      sigma.y = as_labeller(lab_sigma, label_parsed),
      beta = as_labeller(beta_map, label_parsed)
    )
  )+
  theme_bw() +
  scale_color_manual(values = c("black", "grey65")) +
  ylab("Empirical power") +
  xlab(expression(tau)) +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.key.width = unit(1.5, "cm"),
    legend.key.size = unit(4, "line")
  )
