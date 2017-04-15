#################################################################
# Toy example:
#
# p: multivariate normal (0, \Sigma)
# Restrictions: force 0.5 mass of the marginals to lie above some
#               randomly chosen thresholds (chi).
# The above implies:
#                   n = k
#                   f(x) = (x > chi)
#                   m = rep(0.5, n)
#################################################################

library(mcMDISolver)

set.seed(1313) # for reproducibility

n = 9 # number of variables
k = n # number of restrictions
S = 50000 # number of samples

# restrictions
chi = rnorm(n, 1.96, 0.2) # critical values
f = function(x, c = chi){ # vectorized over the size of x (n)
  return(x > c)
}
m = rep(.5, k) # rhs of restrictions

# p-sampler
# p is multivariate normal with constant correlation rho
rho = 0.6 # average correlation
R = matrix(rep(rho, n * n), nrow = n)
diag(R) = rep(1, n)

## SOLVE
start.time = proc.time()
fit_MDI = mcMDISolver::MDI_solve(mvtnorm::rmvnorm(n = S, sigma = R), f = f, m = m)
print(proc.time() - start.time)
print(fit_MDI$x)

## Examine q distribution
# Draw samples using Metropolis algorithm

# log density of q
log_dq = function(x, l = fit_MDI$x){
  return(mvtnorm::dmvnorm(x, sigma = R, log = T) - as.numeric(crossprod(l, c(1, f(x)))))
}

# draw from q
mcmc_S = S # effective number of samples
burn = trunc(S) # number of samples to burn
sample_q = mcmc::metrop(log_dq, initial = chi, nbatch = mcmc_S + burn)
X_q = sample_q$batch[-(1:burn), ]

# histograms
# note the sudden increase in density around chi
par(mfrow = c(3, 3))
for(i in 1:n) hist(X_q[,i], main = paste("chi =", round(chi[i], 2)), xlab = "")
par(mfrow = c(1, 1))

# check restrictions close to zero
rowMeans(apply(X_q, 1, f)) - m
