#' Evaluate target function's Jacobian for equation solver
#'
#' @param L current \eqn{\lambda} estimate of size \eqn{k+1}
#' @param M a list of matrices
#' @param ... currently ignored
#' @return A matrix of size \eqn{(k+1) x (k+1)}
#' @export
J_fun = function(L, M, ...){
  J_list = lapply(M, function(A){
    exp(-sum(L[-1] * A[-1, 1])) * A
    })
  return(exp(-L[1]) * Reduce('+', J_list) / length(J_list))
}

#' Evaluate target function for equation solver
#'
#' @param L current \eqn{\lambda} estimate of size \eqn{k+1}
#' @param M a list of matrices
#' @param m rhs of restrictions
#' @return A vector of size \eqn{k+1}
#' @export
F_fun = function(L, M, m, ...){
  F_list = lapply(M, function(A){
    exp(-sum(L[-1] * A[-1, 1])) * A[, 1]
    })
  return(c(1, m) - exp(-L[1]) * Reduce('+', F_list) / length(F_list))
}

#' Solve MDI with linear equality constraints using a Monte Carlo approach.
#'
#' This function solves a Minimum Discrimination Information problem with linear
#' equality constraints, as described in the package vignette.
#'
#' The non - linear set of equations is solved using package \code{nleqslv}
#' with the Newton method. The rest of the parameters of this solver that
#' are not passed through the \code{control} list are left in their
#' default values.
#'
#' This method assumes that you are able to sample from the \eqn{p} distribution
#' in order to produce the matrix \code{X}. If this is not the case, then
#' a robust alternative is to use an MCMC method (such as the Metropolis algorithm),
#' assuming that you can at least evaluate the density of \eqn{p}.
#'
#' The function \code{f} should take as input a vector \eqn{x} of size \eqn{n} (the
#' dimension of the space where \eqn{x} lives) and return a vector of size
#' \eqn{k} (the number of non-trivial conditions), such that the first component
#' of \eqn{f(x)} corresponds to the value of the first restriction function evaluated
#' at \eqn{x}, the second component to the second restriction, and so on.
#'
#' Keep in mind that if you use the parallelized implementation, and \code{f}
#' depends on some other variables, you need to export those using
#' \code{parallel::clusterExport}.
#'
#' @param X an \eqn{S x n} matrix containing \eqn{S} samples from \eqn{p}
#' @param f a function returning a vector of size \eqn{k} that evaluates
#'  the LHS of restrictions
#' @param m a vector of size \eqn{k} with the RHS of the restrictions
#' @param l_start a guess of the optimal \eqn{k + 1} lambdas. Default is \code{rep(0, k + 1)}
#' @param control a list of parameters passed to \code{nleqslv}
#' @param cl a cluster object for parallelized implementation
#' @return A list as returned by \code{nleqslv}
#' @seealso \code{\link{parallel}}, \code{\link{nleqslv}} for details
#' on this non - linear equation solver, and \code{\link{mcmc}} for an
#' implementation of the Metropolis algorithm.
#' @example /inst/examples/test.R
#' @export
MDI_solve = function(X, f, m, l_start = rep(0, k + 1), control=list(trace = 1), cl = NULL){

  k = length(m) # number of (non-trivial) restrictions

  # check consistency of inputs
  if(length(f(as.vector(X[1, ]))) != k){
    stop(message = "\n\n Arguments X, f and / or m are inconsistent.\n\n")
    return(NULL)
  }

  cat("Creating auxiliary list of matrices M\n")
  M = lapply(split(t(X), rep(1:nrow(X), each = ncol(X))), # convert X to list of rows
             FUN = function(x){
    return(tcrossprod(c(1, f(x))))
  })

  # throw away matrix of samples and take out garbage
  rm(X); gc(verbose = F)

  if(is.null(cl)){
    cat("'cl' argument not present. Using serial implementation\n")
    fit_slv = MDI_solve_ser(M = M, f = f, m = m, l_start = l_start, control = control)
  }else{
    cat("Using parallelized implementation\n")
    fit_slv = MDI_solve_par(M = M, f = f, m = m, l_start = l_start, control = control, cl = cl)
  }

  cat(paste0("\nnleqslv ended with condition ", fit_slv$termcd, ": ", fit_slv$message, "\n\n"))

  return(fit_slv)

}

#' Solve MDI serially.
#'
#' Shouldn't be called directly by user.
#'
#' @inheritParams MDI_solve
#' @return An object of class \code{nleqslv}
MDI_solve_ser = function(M, f, m, l_start, control){

  # solve
  cat("Launching solver\n\n")
  fit_slv = nleqslv::nleqslv(x = l_start, fn = F_fun, jac = J_fun, M = M, m = m,
                             method = "Newton", control = control)

  return(fit_slv)

}
