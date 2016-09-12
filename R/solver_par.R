#' Solve MDI in parallel
#'
#' Shouldn't be called directly by user.
#'
#' @inheritParams MDI_solve
#' @return An object of class \code{nleqslv}
MDI_solve_par = function(X, f, m, l_start = NULL, maxit = 200, cl){

  S = length(X) # number of samples

  # create environment with variables needed for F_fun and J_fun
  cat("Creating auxiliary list of matrices M\n")
  M = parallel::parLapply(cl, X, fun = function(x){
    return(tcrossprod(c(1, f(x))))
  })

  # throw away matrix of samples and take out garbage
  rm(X); gc(verbose = F)

  # export M
  cat("Exporting M to workers\n\n")
  parallel::clusterExport(cl = cl, varlist = list("M"), envir = environment())

  # solve
  cat("Launching solver\n\n")
  fit_slv = nleqslv::nleqslv(x = l_start, fn = F_fun_par, jac = J_fun_par, S = S, m = m,
                             cl = cl, method = "Newton", control=list(maxit = maxit, trace = 1))

  cat(paste0("\nnleqslv ended with condition ", fit_slv$termcd, ": ", fit_slv$message, "\n\n"))

  return(fit_slv)

}

#' Evaluate target function's Jacobian for equation solver in cluster
#'
#' This function shouldn't be called directly by the user.
#'
#' @inheritParams J_fun
#' @param S number of samples
#' @param cl a cluster
#' @return A matrix of size \eqn{(k+1) x (k+1)}
J_fun_par = function(L, S, m, cl){
  tempo = parallel::parLapply(cl, 1:S, function(L, i){
    A = M[[i]]
    exp(-as.numeric(crossprod(L[-1], A[-1, 1]))) * A
    }, L = L)
  return(exp(-L[1]) * Reduce('+', tempo) / S)
}

#' Evaluate target function for equation solver in cluster
#'
#' This function shouldn't be called directly by the user.
#'
#' @inheritParams F_fun
#' @param S number of samples
#' @param cl a cluster
#' @return A vector of size \eqn{k+1}
F_fun_par = function(L, S, m, cl){
  tempo = parallel::parLapply(cl, 1:S, function(L, i){
    A = M[[i]]
    exp(-as.numeric(crossprod(L[-1], A[-1, 1]))) * A[, 1]
    }, L = L)
  return(c(1, m) - exp(-L[1]) * Reduce('+', tempo) / S)
}
