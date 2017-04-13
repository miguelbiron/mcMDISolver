#' Solve MDI with linear equality constraints using a Monte Carlo approach.
#'
#' Same as the other but using estimates of the function and jacobian
#'
#' @param rp a function that takes an integer \code{S} as input and returns \code{S} samples from \eqn{p}
#' @param f a function returning a vector of size \eqn{k} that evaluates
#'  the LHS of restrictions
#' @param m a vector of size \eqn{k} with the RHS of the restrictions
#' @param l_start a guess of the optimal \eqn{k + 1} lambdas. Default is \code{rep(0, k + 1)}
#' @param lr the learning rate
#' @param S the size of the samples to draw from p at each iteration
#' @param eps error tolerance (supremum norm)
#' @param max_iter maximum number of iterations
#' @return A list
#' @example /inst/examples/test.R
#' @export
solve_stoch = function(rp, f, m, l_start = rep(0, k + 1), lr = 0.01, S = 100, eps = 1e-4, max_iter = 1000){

  k = length(m) # number of (non-trivial) restrictions

  # check consistency of inputs
  if(length(f(as.vector(rp(1)))) != k){
    stop(message = "\n\n Arguments rp, f and / or m are inconsistent.\n\n")
    return(NULL)
  }

  # auxiliary variables
  curr_f = rep(2 * eps, k + 1)
  curr_l = l_start
  i = 1

  while(i <= max_iter){

    # check convergence
    err = max(abs(curr_f))
    if(err < eps){
      cat(sprintf("\nConvergence reached after %d iterations", i - 1))
      break
    }

    # print status
    if(i %% 10 == 0){
      cond = sqrt(sum(curr_J ^2)) * sqrt(sum(ginv_J ^2))
      cat(sprintf("\nIteration %d: error = %.4f, cond = %E", i, err, cond))
    }

    # obtain samples and transform into list of matrices
    X = rp(S) # obtain matrix of samples (size S x n)
    M = lapply(split(t(X), rep(1:nrow(X), each = ncol(X))), # convert X to list of rows
               FUN = function(x){
                 return(tcrossprod(c(1, f(x))))
               })

    # take step
    curr_f = F_fun(L = curr_l, M = M, m = m)
    curr_J = J_fun(L = curr_l, M = M)
    ginv_J = MASS::ginv(curr_J)
    h_opt = as.vector(ginv_J %*% ((-1)*curr_f))
    curr_l = curr_l + lr * (h_opt / sqrt(sum(h_opt^2)))
    curr_f = F_fun(curr_l, M, m)
    i = i+1

  }

  return(list(
    lambda = curr_l,
    f_val = curr_f,
    n_iter = i - 1
  ))

}
