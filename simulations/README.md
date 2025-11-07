# Simulation Studies for *Inference on Multiple Quantiles in Regression Models by a Rank-Score Approach*

This folder contains the R scripts used to reproduce the simulation studies presented in **Section 5** of the paper  
**"Inference on Multiple Quantiles in Regression Models by a Rank-Score Approach"**  
by *De Santis, Vesely, and Andreella*.

All simulations rely on the methods implemented in the [`quasar`](https://github.com/angeella/quasar/tree/main/simulations) R package.

---

## Overview

The simulation studies assess:
1. **Type I error control** under the null hypothesis (`β = 0`)
2. **Power behavior** under the alternative hypothesis (`β ≠ 0`)

The scripts generate results corresponding to **Figures 2 and 3** in the paper.

---

## Files

### `error.R`
Simulates data under the **null hypothesis** to assess **Type I error control**.

- Corresponds to **Figure 2 (Section 5)**.  
- Simulates 1000 datasets with `β = 0` under a normal model.  
- Evaluates:
  - **Rank-score tests**
  - **Wald-type tests**
- Reports the empirical rejection rate (proportion of p-values ≤ 0.05).

**Main output:**  
`TypeIcontrol` — a data frame containing the empirical Type I error for each test component.

---

### `power.R`
Simulates data under the **alternative hypothesis** to analyze **power performance**.

- Corresponds to **Figure 3 (Section 5)**.  
- Varies across:
  - Distributions of the response: `t`, `exponential`, `normal`
  - Scale models: `σ_y = 1` or `σ_y = 1 + |X|`
  - Parameter of interest: `β ∈ {0.2, 0.4, 0.6, 0.8, 1.0, 1.2}`
- Computes adjusted p-values using:
  - **Closed testing (rank-score)**
  - **Bonferroni correction**

**Main output:**  
`res` — a data frame reporting the empirical power (proportion of p-values ≤ 0.05) for each configuration.

---

## Requirements

- **R (≥ 4.0)**
- Packages:
  - [`quasar`](https://github.com/angeella/quasar)
  - [`quantreg`](https://cran.r-project.org/package=quantreg)

To install `quasar` from GitHub:
```r
# install.packages("remotes")
remotes::install_github("angeella/quasar")
