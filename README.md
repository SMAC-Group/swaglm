# `swaglm` Overview <img src="man/figures/logo.png" align="right" style="width: 15%; height: 15%"/>

<!-- badges: start  -->
![](https://img.shields.io/github/last-commit/SMAC-Group/swaglm) 
[<img src="https://s-a.github.io/license/img/agpl-3.0.svg" />](https://s-a.github.io/license/?license=agpl-3.0&fullname=Stephan%20Ahlf&year=2015&profile=https://github.com/s-a&projectUrl=https://github.com/s-a/license&projectName=License%20Demo "")
![R-CMD-check](https://github.com/SMAC-Group/swaglm/actions/workflows/R-CMD-check.yaml/badge.svg)
<!-- badges: end -->




## Overview
The `swaglm` package is a fast implementation of the Sparse Wrapper Algorithm (SWAG) for Generalized Linear Models (GLM). SWAG is a meta-learning procedure that combines screening and wrapper methods to efficiently find strong low-dimensional attribute combinations for prediction. Additionally, the package provides a statistical test to assess whether the selected models (learners) extract meaningful information from the data.

## Features
- Efficiently finds a set of low-dimensional learners with high predictive accuracy.
- Follows a forward-step method to iteratively build strong learners.
- Provides a permutation-based statistical test (`swag_test`) to determine if the obtained models capture meaningful structure in the data.
- Uses entropy-based network measures (entropy of frequency and entropy of eigenvalue centrality) to compare SWAG models against randomized models.

Below are instructions on how to install and make use of the `swaglm` package.

## Installation Instructions

The `swaglm` package is currently only available on  GitHub.

``` r
# Install dependencies
install.packages(c("devtools"))

# Install/Update the package from GitHub
devtools::install_github("SMAC-Group/swaglm")

# Install the package with Vignettes/User Guides 
devtools::install_github("SMAC-Group/swaglm", build_vignettes = TRUE)
```


### External `R` libraries

The `swaglm` package relies on a limited number of external libraries, but notably on `Rcpp` and `RcppArmadillo` which require a `C++` compiler for installation, such as for example `gcc`.




## Getting started

```r
library(swaglm)

# Simulated data
n <- 2000
p <- 50
X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = diag(rep(1/p, p)))
beta <- c(-15, -10, 5, 10, 15, rep(0, p-5))
z <- 1 + X %*% beta
pr <- 1 / (1 + exp(-z))
y <- as.factor(rbinom(n, 1, pr))
y <- as.numeric(y) - 1

# Run SWAG
swag_obj <- swaglm(X = X, y = y, p_max = 20, family = binomial(),
                   alpha = 0.15, verbose = TRUE, seed = 123)


# Run statistical test
test_results <- swag_test(swag_obj, B = 50, verbose = TRUE)

# View p-values for both entropy-based measures
print(test_results)
```

## How the statistical test works

The function `swag_test()` performs a permutation test to evaluate whether the selected variables contain meaningful information or are randomly selected.

Null Hypothesis: The selected models are no different from randomly chosen ones.

### Procedure:

- The response variable is shuffled to break its true relationship with predictors.

- SWAG is applied to these shuffled datasets.

- The entropy of variable frequency and eigenvalue centrality is computed for the null models.

- p-values are computed by comparing the SWAG network with these null models.

### Interpretation:

- Small p-value (< 0.05): The selected variables are likely informative.

- Large p-value (≥ 0.05): The selection may be random.


## License

This source code is released under is the GNU AFFERO GENERAL PUBLIC LICENSE (AGPL) v3.0. 

## References

[Molinari, R. et al. SWAG: A Wrapper Method for Sparse Learning (2021) ](https://arxiv.org/abs/2006.12837)


