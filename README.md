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
- Works on top of the `caret` package and follows a forward-step method to iteratively build strong learners.
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


## License

This source code is released under is the GNU AFFERO GENERAL PUBLIC LICENSE (AGPL) v3.0. 

## References

[Molinari, R. et al. SWAG: A Wrapper Method for Sparse Learning (2021) ](https://arxiv.org/abs/2006.12837)


