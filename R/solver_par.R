#' Solve MDI in parallel
#'
#' Shouldn't be called directly by user.
#'
#' @inheritParams MDI_solve
#' @return An object of class \code{nleqslv}
MDI_solve_par = function(M, f, m, l_start, control, cl){

  S = length(M) # number of samples

  # dividing the work among workers
  sched = split(1:S, ceiling(1:S / (S / length(cl))))

  # export M and schedule
  cat("Exporting data to workers\n\n")
  parallel::clusterExport(cl = cl, varlist = list("M", "sched"),
                          envir = environment())

  # solve
  cat("Launching solver\n\n")
  fit_slv = nleqslv::nleqslv(x = l_start, fn = F_fun_par, jac = J_fun_par, S = S, m = m,
                             cl = cl, method = "Newton", control = control)

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
J_fun_par = function(L, S, cl, ...){

  tempo = parallel::parLapply(cl, 1:length(cl), function(i, L){

    # assign a subset of M
    A = M[sched[[i]]] # M should have been exported

    # call the serialized version of the function with a subset of the matrices
    return(J_fun(L, A) * length(A)) # serialized version returns value divided by length

  }, L)

  return(Reduce('+', tempo) / S)

}

#' Evaluate target function for equation solver in cluster
#'
#' This function shouldn't be called directly by the user.
#'
#' @inheritParams F_fun
#' @param S number of samples
#' @param cl a cluster
#' @return A vector of size \eqn{k+1}
F_fun_par = function(L, m, S, cl, ...){
  tempo = parallel::parLapply(cl, 1:length(cl), function(i, L, m){

    # assign a subset of M
    A = M[sched[[i]]] # M should have been exported

    # call the serialized version of the function with a subset of the matrices
    return((F_fun(L, A, m) - c(1, m)) * length(A)) # correct for serial version output

  }, L, m)

  return(c(1, m) + Reduce('+', tempo) / S)
}
