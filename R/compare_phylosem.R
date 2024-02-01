#' Compare phylogenetic structural equation models
#'
#' Fits several phylogenetic structural equation model for further comparison
#'
#' @inheritParams phylosem
#'
#' @param sem_set A named list of structural equation model specifications,
#'        where each element will be passed as argument \code{sem} to
#'        \code{\link{phylosem}}
#' @param ... Additional arguments passed to \code{\link{phylosem}}
#'
#' @return An object (list) of class `compare_phylosem`, containing a list of
#'         output from \code{\link{phylosem}}
#'
#' @export
compare_phylosem <-
function( sem_set,
          tree,
          data,
          family = rep("fixed", ncol(data)),
          covs,
          estimate_ou = FALSE,
          estimate_lambda = FALSE,
          estimate_kappa = FALSE,
          control = phylosem_control(),
          ... ){

  # Initialize object of results
  out = vector( "list", length=length(sem_set) )
  names(out) = names(sem_set)

  # Loop through models
  for( modelI in seq_along(sem_set) ){
    fit = tryCatch( phylosem( sem = sem_set[[modelI]],
    #fit = phylosem( sem = sem_set[[modelI]],
                              tree = tree,
                              data = data,
                              family = family,
                              covs = covs,
                              estimate_ou = estimate_ou,
                              estimate_lambda = estimate_lambda,
                              estimate_kappa = estimate_kappa,
                              control = control,
                              ... ),
                    error = function(e) e )
    if( isTRUE(is(fit, "phylosem")) ){
      out[[modelI]] = fit
    }
  }

  # pass out
  class(out) = "compare_phylosem"
  return(out)
}

#' Extract best fitted model
#'
#' @param x output from \code{\link{compare_phylosem}}
#'
#' @return Returns best model from those fitted using \code{\link{compare_phylosem}}
#'
#' @export
best <- function(x) UseMethod('best')
#' Extract best fitted model
#'
#' @param x output from \code{compare_phylosem}
#'
#' @return Returns best model from those fitted using \code{\link{compare_phylosem}}
#'
#' @method best compare_phylosem
#' @export
best.compare_phylosem <-
function( x ) {
  AICs <- vapply(x, AIC, 1)
  x[[which.min(AICs)]]
}

#' Choose model
#'
#' @param x output from \code{compare_phylosem}
#' @param choice Integer indicating model to extract
#'
#' @return Returns chosen model from those fitted using \code{\link{compare_phylosem}}
#'
#' @export
choice <- function(x, choice) UseMethod('choice')
#' Choose model
#'
#' @param x output from \code{compare_phylosem}
#' @param choice Integer indicating model to extract
#'
#' @return Returns chosen model from those fitted using \code{\link{compare_phylosem}}
#'
#' @method choice compare_phylosem
#' @export
choice.compare_phylosem <-
function( x,
          choice ) {

  x[[choice]]
}

#' Choose model
#'
#' @param x output from \code{compare_phylosem}
#' @param cut_off threshold where any model with delta-AIC greater than this value is excluded from average
#' @param avg_method see \code{\link[phylopath]{average_DAGs}}
#'
#' @return Returns an AIC-weighted average of fitted models from \code{\link{compare_phylosem}} after conversion to format from [phylopath::est_DAG]
#'
#' @export
average <- function(x, cut_off, avg_method) UseMethod('average')
#' Choose model
#'
#' @param x output from \code{compare_phylosem}
#' @param cut_off threshold where any model with delta-AIC greater than this value is excluded from average
#' @param avg_method see \code{\link[phylopath]{average_DAGs}}
#'
#' @return Returns an AIC-weighted average of fitted models from \code{\link{compare_phylosem}} after conversion to format from \code{\link[phylopath]{est_DAG}}
#'
#' @method average compare_phylosem
#' @export
average.compare_phylosem <-
function( x,
          cut_off = 2,
          avg_method = "conditional") {

  AICs <- vapply(x, AIC, 1)
  dAICs <- AICs - min(AICs)

  # calculate Akaike weights
  l <- exp(-0.5 * dAICs)
  w <- l / sum(l)

  selected <- x[dAICs < cut_off]
  #selected_DAGs <- NULL
  #for( i in seq_along(selected) ){
  #  selected_DAGS[[i]] = as( selected[[i]], "fitted_DAG" )
  #}
  selected_DAGs <- lapply(selected, as_fitted_DAG)
  average_DAGs(selected_DAGs, w[dAICs < cut_off], avg_method = avg_method)
}
