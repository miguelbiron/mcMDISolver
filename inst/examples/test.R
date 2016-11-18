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
S = 10000 # number of samples

# restrictions
chi = rnorm(n, 1.96, 0.5) # critical values
f = function(x, c = chi){ # vectorized over the size of x (n)
  return(x > chi)
}
m = rep(.5, k) # rhs of restrictions

# p-sampler
# p is mvnorm with average cor = rho
rho = 0.6 # average correlation
R = diag(n)
R[upper.tri(R)] = pmin(pmax(rnorm(n*(n-1)/2, rho, 0.1), -1), 1)
R[lower.tri(R)] = t(R)[lower.tri(R)]

## SOLVE: serial implementation
start.time = proc.time()
fit_MDI = MDI_solve(mvtnorm::rmvnorm(n = S, sigma = R), f = f, m = m)
print(proc.time() - start.time)

## SOLVE: parallelized implementation

cl = parallel::makeCluster(getOption("cl.cores", 2))

# export additional parameters of f
parallel::clusterExport(cl, list("chi"))

start.time = proc.time()
fit_MDI = MDI_solve(mvtnorm::rmvnorm(n = S, sigma = R), f = f, m = m, cl = cl)
print(proc.time() - start.time)

parallel::stopCluster(cl)

## Examine q distribution
# Draw samples using Metropolis algorithm

# log density of q
log_dq = function(x, l = fit_MDI$x){
  return(mvtnorm::dmvnorm(x, sigma = R, log = T) - as.numeric(crossprod(l, c(1, f(x)))))
}

# draw from q
mcmc_S = 2 * S
burn = trunc(0.1 * mcmc_S) # number of samples to burn
sample_q = mcmc::metrop(log_dq, initial = chi, nbatch = mcmc_S + burn)
X_q = sample_q$batch[-(1:burn), ]

# histograms
# note the sudden increase in density around chi
par(mfrow = c(3, 3))
for(i in 1:n) hist(X_q[,i], main = paste("chi =", round(chi[i], 2)), xlab = "")
par(mfrow = c(1, 1))

# check restrictions close to zero
rowMeans(apply(X_q, 1, f)) - m
