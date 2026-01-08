# Simulation Studies for *Inference on Multiple Quantiles in Regression Models by a Rank-Score Approach*

This folder contains the `R` scripts used to reproduce the simulation studies presented in **Section 5** of the paper  
**“Inference on Multiple Quantiles in Regression Models by a Rank-Score Approach”**  
by *De Santis, Vesely, and Andreella*.

All simulations rely on the methods implemented in the `quasar` `R` package.

---

## Overview

The simulation studies assess:
1. **Type I error control** under the null hypothesis (`β = 0`)
2. **Power behavior** under the alternative hypothesis (`β ≠ 0`)

The scripts generate results corresponding to **Figures 2, 3, 4, and 5** in the paper.

---

## Files

### `error.R`

Simulates data under the **null hypothesis** to assess **Type I error control**.

- Corresponds to **Figure 2 (Section 6)**  
- Generates 1,000 datasets under a normal model with `β = 0`  
- Quantile levels: `{0.1, 0.25, 0.5, 0.75, 0.9}`  
- Compares:
  - **Generalized rank-score tests**
  - **Wald-type tests**
- Reports empirical rejection rates (proportion of p-values ≤ 0.05)

**Main output:**  
`TypeIcontrol` — a data frame containing empirical Type I error rates for each test component.

---

### `power.R`

Simulates data under the **alternative hypothesis** to analyze **power performance**.

- Corresponds to **Figure 3 (Section 6)**  
- Generates 1,000 datasets under a normal model with `β ≠ 0`  
- Quantile levels: `{0.1, 0.25, 0.5, 0.75, 0.9}`  
- Compares **generalized rank-score tests** with:
  - **Base multivariate specification** (`B = A^{-1}`)
  - **Identity weighting** (`B = I_5`)
- Reports empirical power (proportion of p-values ≤ 0.05)

**Main output:**  
`power` — a data frame containing empirical power values for each test component.

---

### `extreme.R`

Simulates data under the **alternative hypothesis** to analyze **power at extreme quantiles**.

- Corresponds to **Figure 4 (Section 6)**  
- Quantile levels: `{0.05, 0.10, 0.15}`  
- Varies:
  - Response distributions: `t`, `skew-normal`, `normal`  
  - Scale specifications: `σ_y = 1`, `σ_y = 1 + |X|`, `σ_y = 1 + |Z|`  
  - Parameter of interest: `β ∈ {0.2, 0.4, 0.6, 0.8, 1.0, 1.2}`  
- Considers three weighting-matrix specifications (dimension `3 × 3`):
  - **Identity weighting** (`B = I_3`)
  - **Inverse diagonal of Δ** (variance-balancing diagonal weighting)
  - **Density-based diagonal weighting**, with diagonal elements proportional to  
    `f(F^{-1}(τ_1)), f(F^{-1}(τ_2)), f(F^{-1}(τ_3))`
- Computes adjusted p-values using:
  - **Closed testing (rank-score)**
  - **Holm–Bonferroni**
  - **Bonferroni**

**Main output:**  
`sim` — a data frame containing, for each simulation configuration:

- **Adjusted p-values for closed testing**
  - `pvC1`, `pvC2`, `pvC3`
- **Holm-adjusted p-values**
  - `pvH1`, `pvH2`, `pvH3`
- **Bonferroni-adjusted p-values**
  - `pvB1`, `pvB2`, `pvB3`

Each triplet corresponds to the three individual hypotheses associated with the quantile levels  
`τ = {0.05, 0.10, 0.15}`.

---

### `decile.R`

Simulates data under the **alternative hypothesis** to analyze **power across decile quantiles**.

- Corresponds to **Figure 5 (Section 5)**  
- Quantile levels: `{0.1, 0.2, …, 0.9}`  
- Varies:
  - Response distributions: `t`, `skew-normal`, `normal`
  - Scale specifications: `σ_y = 1`, `σ_y = 1 + |X|`, `σ_y = 1 + |Z|`
  - Parameter of interest: `β ∈ {0.2, 0.4, 0.6, 0.8, 1.0, 1.2}`
  - Weighting matrices: identity, density-based diagonal, inverse diagonal
- Computes adjusted p-values using:
  - **Closed testing (rank-score)**
  - **Holm–Bonferroni**
  - **Bonferroni**

**Main output:**  
`sim` — a data frame containing, for each simulation configuration:

- **Closed testing adjusted p-values**
  - `pvC1`, …, `pvC9`
- **Holm-adjusted p-values**
  - `pvH1`, …, `pvH9`
- **Bonferroni-adjusted p-values**
  - `pvB1`, …, `pvB9`


Each triplet (`*1`, …, `*9`) corresponds to the nine individual hypotheses associated with  
`τ = {0.1, 0.2, …, 0.9}`.
---

## Requirements

- **R (≥ 4.0)**
- Packages:
  - `quasar`
  - `quantreg`
