<img src="sticker.svg" align="right" alt="" width="200" />

# quasar

`quasar` is the package developed to have valid inference when multiple quantile regressions are fitted. 

# Toy example 

`remotes::install_github("angeella/quasar")
`library(quasar)`

`p <- 3; Sigma <- diag(p)`

`dat_n <- simulateData(200, beta=0, gamma=c(0.2,-0.1), mu=4, Sigma=Sigma, sigma.y=0.5, distribution="t")`

`mod <- rq(y ~ X + Z1 + Z2, tau = c(0.1, 0.25, 0.5, 0.75, 0.9), data = dat_n)`

`closedTesting(mod, X = "X")`
