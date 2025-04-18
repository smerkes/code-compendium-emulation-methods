# CODE COMPENDIUM: MULTIVARIATE AND FUNCTIONAL OUTPUT EMULATION

This repository contains R code for emulating computer model outputs with multivariate and functional structure using Gaussian Process models. It includes three main source files supporting parameter estimation, efficient posterior computation via Kronecker product structures, and simulation from conditional distributions. An accompanying R Markdown file and revised PDF chapter will guide users through the code and underlying methodology.

## PDF Chapter and Rmd files related 

  * **mvEmulatorChapterCode.pdf**: A follow-along for Chapter 7 in Adventures in Space-Time Modeling with Maike and Dave.
  * **mvEmulatorChapterCode.Rmd**: *in development* to allow the audience to follow along more easily. 

## .R Function Files
Below is a list of the individual functions within each R script to help guide users in understanding where specific functions originate.

`Cov_functions.R` contains the following functions: 

  * `eudis`: calculate Euclidean distance
  * `get_distmat`: construct a distance matrix (distance is Euclidean distance)
  * `rmultnorm`: generates realizations from a multivariate normal (MVN) distribution
  * `rmultnorm_SVD`: efficiently generates realizations from an MVN normal distribution where the matrix is represented as `kronecker(cov1, cov2)`

`Estimation_1D.R` contains the following functions and library:

  * Uses `dfoptim` library in R: (https://cran.r-project.org/web/packages/dfoptim/index.html)
  * `Fast_kronecker_quadratic`: efficiently carries out the quadratic equation for MVN likelihood as well as the log determinant
  * `LogLik_MVN`: this function calls Fast_kronecker_quadratic() and calculates the log likelihood
  * `Estimate_params`: estimate parameters for 1D estimation of marginal variance, scales, and nugget given data z. This function assumes the udnerlying distribution is MVN(0, `kronecker(C, R)`)

`Estimation_Basis.R` contains the following functions:

  * `ML_pcEMU`: function that estimates $\phi_j$'s and $\sigma^2_j$'s using maximum likelihood
  * `Estimate_params_pcEm`: main function for optimization - called by `ML_pcEMU()`
  * `pre_process`: Pre-process computer model output to facilitate parameter estimation. This function is called by `ML_pcEMU()`
  * `get_wstar_distr_preds`: estimate phi and $\sigma^2$ for $w^*_j$
