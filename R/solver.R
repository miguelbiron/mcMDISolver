#' Evaluate target function's Jacobian for equation solver serially
#'
#' This function shouln't be called directly by the user.
#'
#' @param L current \eqn{\lambda} estimate of size \eqn{k+1}
#' @return A matrix of size \eqn{(k+1) x (k+1)}
J_fun = function(L, M, m){
  tempo = lapply(M, function(A){exp(-as.numeric(crossprod(L[-1], A[-1, 1]))) * A})
  return(exp(-L[1]) * Reduce('+', tempo) / length(tempo))
}

#' Evaluate target function for equation solver serially
#'
#' This function shouln't be called directly by the user.
#'
#' @param L current \eqn{\lambda} estimate of size \eqn{k+1}
#' @return A vector of size \eqn{k+1}
F_fun = function(L, M, m){
  tempo = lapply(M, function(A){exp(-as.numeric(crossprod(L[-1], A[-1, 1]))) * A[, 1]})
  return(c(1, m) - exp(-L[1]) * Reduce('+', tempo) / length(tempo))
}

#' Solve MDI with linear constraints using a Monte Carlo approach.
#'
#' Solve MDI with linear constraints using a Monte Carlo approach.
#'
#' @param X an \eqn{S x n} matrix containing samples from \eqn{p}
#' @param f a function returning a vector of size \eqn{k} that evaluates the LHS of restrictions
#' @param m a vector of size \eqn{k} with the RHS of the restrictions
#' @param l_start a guess of the optimal \eqn{k + 1} lambdas. Default is \eqn{0_{k+1}}
#' @param maxit maximum number of iterations for \code{nleqslv}. Default is 200
#' @param cl a cluster object for parallelized implementation.
#' @return An object of class \code{nleqslv}
#' @export
MDI_solve = function(X, f, m, l_start = NULL, maxit = 200, cl = NULL){

  # check consistency of inputs
  if(length(f(as.vector(X[1, ]))) != length(m)){
    stop(message = "\n\n Arguments X, f and / or m are inconsistent.\n\n")
    return(NULL)
  }

  k = length(m) # number of (non-trivial) restrictions

  # assign a default l_start
  if(is.null(l_start)){
    l_start = rep(0, k + 1)
  }

  # convert X to list of rows
  X = split(t(X), rep(1:nrow(X), each = ncol(X)))
  gc(verbose = F)

  if(is.null(cl)){
    cat("'cl' argument not present. Using serial implementation\n")
    return(MDI_solve_ser(X = X, f = f, m = m, l_start = l_start, maxit = 200))
  }else{
    cat("Using parallelized implementation\n")
    return(MDI_solve_par(X = X, f = f, m = m, l_start = l_start, maxit = 200, cl = cl))
  }

}

#' Solve MDI serially.
#'
#' Shouldn't be called directly by user.
#'
#' @inheritParams MDI_solve
#' @return An object of class \code{nleqslv}
MDI_solve_ser = function(X, f, m, l_start = NULL, maxit = 200){

  # create environment with variables needed for F_fun and J_fun
  cat("Creating auxiliary list of matrices M\n")
  M = lapply(X, FUN = function(x){
    return(tcrossprod(c(1, f(x))))
  })

  # throw away matrix of samples and take out garbage
  rm(X); gc(verbose = F)

  # solve
  cat("Launching solver\n\n")
  fit_slv = nleqslv::nleqslv(x = l_start, fn = F_fun, jac = J_fun, M = M, m = m,
                             method = "Newton", control=list(maxit = maxit, trace = 1))

  cat(paste0("\nnleqslv ended with condition ", fit_slv$termcd, ": ", fit_slv$message, "\n\n"))

  return(fit_slv)

}
