# iivDecomp: Decompose Intraindividual Variability

[![R](https://img.shields.io/badge/R-%3E%3D%204.0-blue)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

An R package for decomposing observed intraindividual variability (IIV) in intensive longitudinal data into four components: true construct-related variability, measurement error, careless responding (insufficient effort responding; IER), and context effects.

## Installation

```r
# Install from GitHub
remotes::install_github("user/iivDecomp")
```

## Quick Start

```r
library(iivDecomp)

# Generate simulated data with 15% IER contamination
dat <- iiv_generate(
  N = 200, T = 30, K = 8,
  sigma_tau = 0.3, phi = 0.5,
  p_IER = 0.15, IER_type = "random",
  seed = 42
)

# Run the decomposition
res <- iiv_decompose(dat)
print(res)

# Compare all five methods
iiv_compare(dat)
```

## Core Functions

| Function | Description |
|----------|-------------|
| `iiv_decompose()` | EM-based variance decomposition (M4 method) |
| `iiv_fit_lmer()` | Multilevel model for within-person variance |
| `iiv_compare()` | Compare M1-M4 methods on one dataset |
| `iiv_generate()` | Simulate intensive longitudinal data |
| `iiv_detect_longstring()` | Longstring IER detection |
| `iiv_detect_mahalanobis()` | Mahalanobis distance IER detection |
| `iiv_aggregate()` | Classical aggregation (M1) |

## Methods

The package implements four IIV estimation approaches:

- **M1 (Classical Aggregation):** Treats all within-person variance as noise
- **M2 (Multilevel/DSEM):** Treats all within-person residual variance as signal
- **M3a (Longstring Screening):** Binary IER detection + refit
- **M3b (Mahalanobis Screening):** Outlier detection + refit
- **M4 (EM Decomposition):** Continuous variance decomposition via EM mixture model

## Key Findings

From a simulation study (5,040 conditions, N = 50-200, T = 10-100, K = 4-8):

- **M4 reduces RMSE by 50% overall** vs. M2, and up to 80% in high-IER conditions
- At 30% IER contamination, **M4 bias is only 1/5 of M2 bias**
- **No single binary IER detector** works for all IER types (random vs. fixed)
- **Robust to missing data** (MCAR/MAR at 15-30%) and model misspecification

## Data Format

Input data should have columns: `id`, `time`, `item`, `y` (observed response). Optionally include `ys` (continuous pre-discretization values) for improved measurement error estimation.

```r
head(dat)
#   id time item y       ys tau_true   ier missing
# 1  1    1    1 3 3.621543   0.1234 FALSE   FALSE
# 2  1    1    2 4 3.891234   0.1234 FALSE   FALSE
# ...
```

## Reference



## License

MIT
