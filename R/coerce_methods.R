#' @title Convert phylosem to phylopath output
#'
#' @description Convert output from package phylosem to phylopath
#'
#' @param object Output from \code{\link{phylosem}}
#'
#' @return Convert output to format supplied by \code{\link[phylopath]{est_DAG}}
#'
#' @export
as_fitted_DAG <-
function( object ){

  # extract and name identical to output from est_DAG
  out = list(
    coef = t(object$report$Rho_jj),
    se = t(as.list(object$opt$SD, what="Std. Error", report=TRUE)$Rho_jj)
  )
  dimnames(out$coef) = dimnames(out$se) = list( colnames(object$data), colnames(object$data) )

  # pass out
  class(out) = "fitted_DAG"
  return(out)
}

#' @title Convert phylosem to sem output
#'
#' @description Convert output from package phylosem to output from package sem
#'
#' @param object Output from \code{\link{phylosem}}
#'
#' @return Output converted to format supplied by \code{\link[sem]{sem}}
#'
#' @export
as_sem <-
function( object ){

  Sprime = object$report$V_jj
    rownames(Sprime) = colnames(Sprime) = colnames(object$data)
  out = sem( object$SEM_model,
             S = Sprime,
             N = nrow(object$data) )

  # pass out
  return(out)
}

#' @title Convert phylosem to phylo4d
#'
#' @description Convert output from package phylosem to phylo4d object from package phylobase
#'
#' @details
#' This package is intended to for use in using plots assocaited with package sem,
#' e.g., using package plotSEM \code{semPlot::semPlotModel}
#'
#' @param object Output from \code{\link{phylosem}}
#' @param what Select what to convert (Estimate / Std. Error).
#'
#' @return phylosem output to converted format supplied by \code{\link[phylobase]{phylo4d}}
#'
#' @export
as_phylo4d <-
function( object,
          what = c("Estimate", "Std. Error") ){

  #
  what = match.arg(what)
  traits = as.list( object$opt$SD, what=what )$x_vj
  colnames(traits) = colnames(object$data)
  out = phylo4d( x=object$tree, all.data=traits )

  # pass out
  return(out)
}
