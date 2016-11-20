
<!-- README.md is generated from README.Rmd. Please edit that file -->
mcMDISolver
===========

This package solves MDI estimation problems with linear equality constraints, using a Monte Carlo approach.

Installation
------------

1.  Install package `devtools`
2.  Run `devtools::install_github("miguelbiron/mcMDISolver")`

Usage
-----

The only function that you will use is `MDI_solve`, feeding it with

-   A matrix `X` with samples from your prior distribution *p*.
-   A vectorized function `f` that evaluates all equality constraints simultaneously.
-   A vector `m` with the right hand side of the constraints.

The output corresponds to a list as returned by `nleqslv::nleqslv`, the non-linear equation solver used.

Toy example
-----------

In this toy example we will shift the mean of a 9-dimensional multivariate normal distribution.

``` r
library(mcMDISolver)
library(mvtnorm) # provides multivariate normal utils
library(mcmc) # to draw samples from the posterior

set.seed(1313) # for reproducibility

n = 9 # number of variables
k = n # number of restrictions equal to the number of variables
S = 50000 # number of samples to draw from prior distribution

# restrictions function
# input: vector of size n. output: vector of size k
f = function(x){ 
  return(x) # first moment
}
m = rep(1, k) # rhs of restrictions. shift mean from 0 to 1

# original mean vector
mu_0 = rnorm(n)
print(mu_0) # inspect the mean
#> [1]  1.1893558  0.5513848  0.3536343 -0.7039707  0.1184256 -0.4236903
#> [7] -0.2177472  0.5442795  1.9105651
```

Now we use the solver function to find the Lagrange multipliers

``` r
## SOLVE
fit_MDI = MDI_solve(
  X = rmvnorm(n = S, mean = mu_0),
  f = f,
  m = m)
#> Creating auxiliary list of matrices M
#> 'cl' argument not present. Using serial implementation
#> Launching solver
#> 
#>   Algorithm parameters
#>   --------------------
#>   Method: Newton  Global strategy: double dogleg (initial trust region = -2)
#>   Maximum stepsize = 1.79769e+308
#>   Scaling: fixed
#>   ftol = 1e-08 xtol = 1e-08 btol = 0.001 cndtol = 1e-12
#> 
#>   Iteration report
#>   ----------------
#>   Iter         Jac     Lambda      Eta     Dlt0     Dltn         Fnorm   Largest |f|
#>      0                                                    4.421580e+00  1.696858e+00
#>      1  N(9.7e-03) N            0.3443   4.3408   0.4341  2.706918e+04  7.732667e+01
#>      1             W   0.2652   0.3443   0.4341   0.8682  3.247539e+00  1.391038e+00
#>      2  N(1.2e-02) W   0.5779   0.3054   0.8682   1.7363  1.790775e+00  1.064393e+00
#>      3  N(1.5e-02) P            0.2640   1.7363   0.1736  4.633756e+01  3.367868e+00
#>      3             W   0.1605   0.2640   0.1736   0.3473  1.483226e+00  9.813941e-01
#>      4  N(1.4e-02) W   0.7964   0.2614   0.3473   0.6945* 9.159142e-01  7.686050e-01
#>      4             P            0.2614   0.6945   1.3890  6.436121e-01  7.555778e-01
#>      5  N(1.3e-02) N            0.4045   1.3246   1.3246  5.602304e-01  5.257754e-01
#>      6  N(1.1e-02) N            0.9583   0.9187   1.8373  2.325957e-02  1.178846e-01
#>      7  N(9.0e-03) N            0.8957   0.4401   0.8802  3.584287e-04  1.442137e-02
#>      8  N(9.2e-03) N            0.9043   0.0632   0.1264  1.216307e-07  2.654283e-04
#>      9  N(8.9e-03) N            0.9119   0.0011   0.0022  1.248698e-14  8.657665e-08
#>     10  N(8.9e-03) N            0.9098   0.0000   0.0000  2.082901e-28  1.088019e-14
#> 
#> nleqslv ended with condition 1: Function criterion near zero
```

Let us take a look at the resulting vector of lagrange multipliers

``` r
print(fit_MDI$x)
#>  [1]  1.88397754  0.04245289 -0.59594223 -0.66858487 -1.90745322
#>  [6] -0.92765604 -1.51250213 -1.24180568 -0.47187529  0.89635811
```

Finally, let us check to see that we have actually shifted the mean to where we wanted. To do this, we will sample values from the posterior distribution using the Metropolis algorithm.

``` r
# Draw samples using Metropolis algorithm
# log density of q
log_dq = function(x, l = fit_MDI$x){
  return(dmvnorm(x, mean = mu_0, log = T) - as.numeric(crossprod(l, c(1, f(x)))))
}

# draw from q
mcmc_S = 50000
burn = trunc(0.1 * mcmc_S) # number of samples to burn
sample_q = metrop(log_dq, initial = m, nbatch = mcmc_S + burn) # samples should be around te new mean -> m
X_q = sample_q$batch[-(1:burn), ]

# check means equal to 1
print(colMeans(X_q))
#> [1] 1.156222 1.153157 1.032946 1.220801 1.036348 1.063682 0.986597 1.022655
#> [9] 1.033310
```

The means are close to the desired output. Increase the size of the sample passed to `MDI_solve` for better accuracy,

Additional information
----------------------

Run `?MDI_solve` for a different example. Also, check the package vignette to get a better understanding of the method used.
