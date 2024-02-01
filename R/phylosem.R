#' @title Fit phylogenetic structural equation model
#'
#' @description Fits a phylogenetic structural equation model
#'
#' @inheritParams sem::specifyModel
#'
#' @param sem structural equation model structure, passed to either \code{\link[sem]{specifyModel}}
#'        or \code{\link[sem]{specifyEquations}} and then parsed to control
#'        the set of path coefficients and variance-covariance parameters
#' @param tree phylogenetic structure, using class \code{\link[ape]{as.phylo}}
#' @param data data-frame providing variables being modeled.  Missing values are inputted
#'        as NA.  If an SEM includes a latent variable (i.e., variable with no available measurements)
#'        then it still must be inputted as a column of \code{data} with entirely NA values.
#' @param family Character-vector listing the distribution used for each column of \code{data}, where
#'        each element must be \code{fixed}, \code{normal}, \code{binomial}, or \code{poisson}.
#'        \code{family="fixed"} is default behavior and assumes that a given variable is measured exactly.
#'        Other options correspond to different specifications of measurement error.
#' @param estimate_ou Boolean indicating whether to estimate an autoregressive (Ornstein-Uhlenbeck)
#'        process using additional parameter \code{lnalpha},
#'        corresponding to the \code{model="OUrandomRoot"} parameterization from \pkg{phylolm}
#'        as listed in \doi{10.1093/sysbio/syu005}
#' @param estimate_lambda Boolean indicating whether to estimate additional branch lengths for
#'        phylogenetic tips (a.k.a. the Pagel-lambda term) using additional parameter \code{logitlambda}
#' @param estimate_kappa Boolean indicating whether to estimate a nonlinear scaling of branch
#'        lengths (a.k.a. the Pagel-kappa term) using additional parameter \code{lnkappa}
#' @param data_labels For each row of \code{data}, listing the corresponding name from
#'        \code{tree$tip.label}.  Default pulls \code{data_labels} from \code{rownames(data)}
#' @param tmb_inputs optional tagged list that overrides the default constructor
#'        for TMB inputs (use at your own risk)
#' @param control Output from \code{\link{phylosem_control}}, used to define user
#'        settings, and see documentation for that function for details.
#'
#' @details
#' Note that parameters \code{logitlambda}, \code{lnkappa}, and \code{lnalpha} if estimated are each estimated as having a single value
#'      that applies to all modeled variables.
#'      This differs from default behavior in \pkg{phylolm}, where these parameters only apply to the "response" and not "predictor" variables.
#'      This also differs from default behavior in \pkg{phylopath}, where a different value is estimated
#'      in each call to \pkg{phylolm} during the d-separation estimate of path coefficients. However, it is
#'      consistent with default behavior in \pkg{Rphylopars}, and estimates should be comparable in that case.
#'      These additional parameters are estimated with unbounded support, which differs somewhat from default
#'      bounded estimates in \pkg{phylolm}, although parameters should match if overriding \pkg{phylolm} defaults
#'      to use unbounded support.  Finally, \code{phylosem} allows these three parameters to be estimated in any
#'      combination, which is expanded functionality relative to the single-option functionality in \pkg{phylolm}.
#'
#' Also note that \pkg{phylopath} by default uses standardized coefficients.  To achieve matching parameter estimates between
#'      \pkg{phylosem} and \pkg{phylopath}, standardize each variable to have a standard deviation of 1.0 prior to fitting with \pkg{phylosem}.
#'
#' @importFrom stats AIC na.omit nlminb optimHess plogis pnorm rnorm
#' @importFrom sem sem pathDiagram specifyModel specifyEquations
#' @importFrom phylopath average_DAGs coef_plot
#' @importFrom phylobase phylo4d
#' @importFrom ape Ntip node.depth.edgelength rtree
#' @importFrom TMB compile dynlib MakeADFun sdreport
#' @importFrom methods is
#'
#' @return
#' An object (list) of class `phylosem`. Elements include:
#' \describe{
#' \item{data}{Copy of argument \code{data}}
#' \item{SEM_model}{SEM model parsed from \code{sem} using \code{\link[sem]{specifyModel}} or \code{\link[sem]{specifyEquations}}}
#' \item{obj}{TMB object from \code{\link[TMB]{MakeADFun}}}
#' \item{tree}{Copy of argument \code{tree}}
#' \item{tmb_inputs}{The list of inputs passed to \code{\link[TMB]{MakeADFun}}}
#' \item{opt}{The output from \code{\link[stats]{nlminb}}}
#' \item{sdrep}{The output from \code{\link[TMB]{sdreport}}}
#' \item{report}{The output from \code{obj$report()}}
#' \item{parhat}{The output from \code{obj$env$parList()} containing maximum likelihood estimates and empirical Bayes predictions}
#' }
#'
#' @references
#' **Introducing the package, its features, and comparison with other software
#' (to cite when using phylosem):**
#'
#' Thorson, J. T., & van der Bijl, W. (In press). phylosem: A fast and simple
#' R package for phylogenetic inference and trait imputation using phylogenetic
#' structural equation models. Journal of Evolutionary Biology.
#' \doi{10.1111/jeb.14234}
#'
#' *Statistical methods for phylogenetic structural equation models*
#'
#' Thorson, J. T., Maureaud, A. A., Frelat, R., Merigot, B., Bigman, J. S., Friedman,
#' S. T., Palomares, M. L. D., Pinsky, M. L., Price, S. A., & Wainwright, P. (2023).
#' Identifying direct and indirect associations among traits by merging phylogenetic
#' comparative methods and structural equation models. Methods in Ecology and Evolution,
#' 14(5), 1259-1275. \doi{10.1111/2041-210X.14076}
#'
#' *Earlier development of computational methods, originally used for phlogenetic factor analysis:*
#'
#' Thorson, J. T. (2020). Predicting recruitment density dependence and intrinsic growth rate for all fishes
#' worldwide using a data-integrated life-history model. Fish and Fisheries, 21(2),
#' 237-251. \doi{10.1111/faf.12427}
#'
#' Thorson, J. T., Munch, S. B., Cope, J. M., & Gao, J. (2017). Predicting life
#' history parameters for all fishes worldwide. Ecological Applications, 27(8),
#' 2262-2276. \doi{10.1002/eap.1606}
#'
#' *Earlier development of phylogenetic path analysis:*
#'
#' van der Bijl, W. (2018). phylopath: Easy phylogenetic path analysis in
#' R. PeerJ, 6, e4718. \doi{10.7717/peerj.4718}
#'
#' von Hardenberg, A., & Gonzalez-Voyer, A. (2013). Disentangling
#' evolutionary cause-effect relationships with phylogenetic confirmatory
#' path analysis. Evolution; International Journal of Organic Evolution,
#' 67(2), 378-387. \doi{10.1111/j.1558-5646.2012.01790.x}
#'
#' *Interface involving SEM `arrow notation` is repurposed from:*
#'
#' Fox, J., Nie, Z., & Byrnes, J. (2020). Sem: Structural equation models.
#' R package version 3.1-11. \url{https://CRAN.R-project.org/package=sem}
#'
#' *Coercing output to phylo4d depends upon:*
#'
#' Bolker, B., Butler, M., Cowan, P., de Vienne, D., Eddelbuettel, D., Holder, M.,
#' Jombart, T., Kembel, S., Michonneau, F., & Orme, B. (2015). phylobase:
#' Base package for phylogenetic structures and comparative data. R Package Version 0.8.0.
#' \url{https://CRAN.R-project.org/package=phylobase}
#'
#' *Laplace approximation for parameter estimation depends upon:*
#'
#' Kristensen, K., Nielsen, A., Berg, C. W., Skaug, H., & Bell, B. M. (2016).
#' TMB: Automatic differentiation and Laplace approximation. Journal of Statistical Software,
#' 70(5), 1-21. \doi{10.18637/jss.v070.i05}
#'
#' @examples
#' # Load data set
#' data(rhino, rhino_tree, package="phylopath")
#'
#' # Run phylosem
#' model = "
#'   DD -> RS, p1
#'   BM -> LS, p2
#'   BM -> NL, p3
#'   NL -> DD, p4
#' "
#' psem = phylosem( sem = model,
#'           data = rhino[,c("BM","NL","DD","RS","LS")],
#'           tree = rhino_tree )
#'
#' # Convert and plot using phylopath
#' library(phylopath)
#' my_fitted_DAG = as_fitted_DAG(psem)
#' coef_plot( my_fitted_DAG )
#' plot( my_fitted_DAG )
#'
#' # Convert to phylo4d to extract estimated traits and Standard errors
#' # for all ancestors and tips in the tree.
#' # In this rhino example, note that species are labeled s1-s100
#' # and ancestral nodes are not named.
#' (traits_est = as_phylo4d(psem))
#' (traits_SE = as_phylo4d(psem, what="Std. Error"))
#'
#' # Convert to sem and plot
#' library(sem)
#' my_sem = as_sem(psem)
#' pathDiagram( model = my_sem,
#'                   style = "traditional",
#'                   edge.labels = "values" )
#' effects( my_sem )
#'
#' # Plot using semPlot
#' if( require(semPlot) ){
#'   myplot = semPlotModel( my_sem )
#'   semPaths( my_sem,
#'                    nodeLabels = myplot@Vars$name )
#' }
#'
#' @useDynLib phylosem, .registration = TRUE
#' @export
phylosem <-
function( sem,
          tree,
          data,
          family = rep("fixed", ncol(data)),
          covs = colnames(data),
          estimate_ou = FALSE,
          estimate_lambda = FALSE,
          estimate_kappa = FALSE,
          data_labels = rownames(data),
          tmb_inputs = NULL,
          control = phylosem_control() ){

  # Function that converts SEM model to a RAM, see `?sem` for more context
  build_ram = function( model, vars ){
    vars = sapply( vars, FUN=function(char){gsub("-", "", gsub(" ", "", char))} )
    n.paths = nrow(model)
    par.names = model[, 2]
    startvalues = model[,3]

    # EXCERPT FROM `getAnywhere("sem.semmod")`
    heads = from = to = rep(0, n.paths)
    for (p in 1:n.paths) {
      #path = sem:::parse.path(model[p, 1])
      path = parse_path(model[p, 1])
      heads[p] = abs(path$direction)
      to[p] = path$second
      from[p] = path$first
      if (path$direction == -1) {
        to[p] = path$first
        from[p] = path$second
      }
    }
    missing_vars = setdiff( c(from,to), vars )
    if( length(missing_vars) > 0 ) stop( "Check `build_ram`:", paste0(missing_vars,sep=", ") )

    ram = data.frame(matrix(0, nrow=p, ncol=5))
    pars = na.omit(unique(par.names))
    ram[, 1] = heads
    ram[, 2] = apply(outer(vars, to, "=="), 2, which)
    ram[, 3] = apply(outer(vars, from, "=="), 2, which)
    par.nos = apply(outer(pars, par.names, "=="), 2, which)
    if(length(par.nos) > 0){
      ram[, 4] = unlist(lapply(par.nos, function(x) if (length(x)==0){0}else{x}))
    }
    ram[, 5] = startvalues
    colnames(ram) = c("heads", "to", "from", "parameter", "start")
    return(ram)
  }

  # Errors / warnings;  not using `assertthat::assert_that` to avoid dependency
  #if( estimate_ou==FALSE & fixed_root==TRUE ){
  #  stop("`fixed_root=TRUE` is only applicable when `estimate_ou=TRUE`")
  #}
  #if( measurement_error==TRUE & estimate_lambda==TRUE ){
  #  stop("measurement errors using `measurement_error=TRUE` are confounded with `estimate_lambda=TRUE`")
  #}
  # General error checks
  if( isFALSE(is(control, "phylosem_control")) ) stop("`control` must be made by `phylosem_control()`")

  if( isFALSE(is(tree, "phylo")) ){
    stop("Check `tree` input")
  }
  if( !("edge.length" %in% names(tree)) ){
    stop("`tree` must include `edge.length` slot")
  }
  familycode_j = sapply( tolower(family), FUN=switch, "fixed"=0, "normal"=1, "norm"=1, "binomial"=2, "binom"=2, "poisson"=3, "pois"=3, "gamma"=4, NA )
  if( any(is.na(familycode_j)) ) stop("Check `family`")
  if( any(is.nan(as.matrix(data))) ) stop("Please remove `NaN` values from data, presumably switching to `NA` values")

  #
  SEM_model = tryCatch(
    specifyModel( text=sem, exog.variances=TRUE, endog.variances=TRUE, covs=covs, quiet=control$quiet ),
    error = function(e) e
  )
  if( isFALSE(is(SEM_model,"semmod")) ){
    SEM_model = tryCatch(
      specifyEquations( text=sem, exog.variances=TRUE, endog.variances=TRUE, covs=covs ),
      error = function(e) e
    )
  }
  if( isFALSE(is(SEM_model,"semmod")) ){
    stop("Must supply either input for `sem::specifyModel` or `sem::specifyEquations`")
  }
  RAM = build_ram( SEM_model, colnames(data) )

  #
  n_tip = Ntip(tree)
  vroot = n_tip + 1
  edge_ez = tree$edge  # parent, child
  length_e = tree$edge.length
  height_v = node.depth.edgelength(tree)
  # seems like length_e + height_v should be equal in an ultrametric tree
  if(vroot %in% edge_ez[,2]) stop("Check for problems")
  if(any(length_e==0)) stop("`tree` contains an edge with length of zero; please fix")

  # associate each datum with tree
  v_i = match( data_labels, c(tree$tip.label,tree$node.label) )
  if(any(is.na(v_i))){
    stop("Check that all `data_labels` are present in `c(tree$tip.label,tree$node.label)`")
  }

  #
  if( is.null(tmb_inputs) ){
    # Build data
    data_list = list( "n_tip" = n_tip,
                      "edge_ez" = edge_ez - 1,
                      "length_e" = length_e,
                      "RAM" = as.matrix(RAM[,1:4]),
                      "RAMstart" = as.numeric(RAM[,5]),
                      "estimate_ou" = estimate_ou,
                      "estimate_lambda" = estimate_lambda,
                      "estimate_kappa" = estimate_kappa,
                      "height_v" = height_v,
                      "y_ij" = as.matrix(data),
                      "v_i" = v_i - 1,
                      "familycode_j" = familycode_j )

    # Build parameters
    rmatrix = function(nrow, ncol) matrix(rnorm(nrow*ncol),nrow=nrow,ncol=ncol)
    parameters_list = list( "beta_z" = rep(0.1, max(RAM[,4])),
                            "lnsigma_j" = rep(0,ncol(data)),
                            "lnalpha" = log(1),
                            "logitlambda" = plogis(2),
                            "lnkappa" = log(1),
                            "x_vj" = 0.1 * rmatrix( nrow=nrow(edge_ez)+1, ncol=ncol(data) ),
                            "xbar_j" = rep(0,ncol(data)) )
    # Build map
    map_list = list()
    # Start off map_list$x_vj, which has multiple constraints
    map_list$x_vj = array( 1:prod(dim(parameters_list$x_vj)), dim=dim(parameters_list$x_vj) )
    # Turn off root for any variable with no measurements (i.e., latent-variables are defined with fixed mean)
    map_list$x_vj[data_list$n_tip+1,] = ifelse( colSums(!is.na(data))==0, NA, map_list$x_vj[data_list$n_tip+1,] )

    # Settings
    map_list$lnsigma_j = 1:length(parameters_list$lnsigma_j)
    for( j in 1:ncol(data)){
      # Turn off SD of measurement error
      if( familycode_j[j] %in% c(0,2,3) ){
        map_list$lnsigma_j[j] = NA
      }
      if( familycode_j[j] == 0 ){
        # Fix random-effects for tips at their observed values
        parameters_list$x_vj[v_i,j] = ifelse( is.na(data_list$y_ij[,j]), parameters_list$x_vj[v_i,j], data_list$y_ij[,j] )
        map_list$x_vj[v_i,j] = ifelse( is.na(data_list$y_ij[,j]), map_list$x_vj[v_i,j], NA )
      }
    }
    if( estimate_ou==FALSE ){
      map_list$lnalpha = factor(NA)
    }
    if( estimate_ou==FALSE ){
      map_list$xbar_j = factor(rep( NA, length(parameters_list$xbar_j) ))
    }else{
      map_list$xbar_j = factor( ifelse(colSums(!is.na(data))==0, NA, 1:ncol(data)) )
    }
    if( estimate_lambda==FALSE ){
      map_list$logitlambda = factor(NA)
    }
    if( estimate_kappa==FALSE ){
      map_list$lnkappa = factor(NA)
    }

    # wrap up map_list$x_vj
    map_list$x_vj = factor(map_list$x_vj)
    map_list$lnsigma_j = factor(map_list$lnsigma_j)

    # Build random
    random = c("x_vj")

    # Bundle
    tmb_inputs = list( map_list=map_list, parameters_list=parameters_list, data_list=data_list, random=random )
  }else{
    if(!all(c() %in% names(tmb_inputs)) ){
      stop("Check contents of `tmb_inputs`")
    }
  }

  # Build TMB object
  obj = MakeADFun( data = tmb_inputs$data_list,
                        parameters = tmb_inputs$parameters_list,
                        map = tmb_inputs$map_list,
                        random = tmb_inputs$random,
                        DLL = "phylosem" )
  if(control$quiet==FALSE) list_parameters(obj)
  results = list( "data" = data,
                  "SEM_model" = SEM_model,
                  "obj" = obj,
                  "call" = match.call(),
                  "tree" = tree,
                  "tmb_inputs" = tmb_inputs )

  # Export stuff
  if( control$run_model==FALSE ){
    return( results )
  }

  #
  obj$env$beSilent()       # if(!is.null(Random))
  #results$opt = fit_tmb( obj,
  #                        quiet = quiet,
  #                        control = list(eval.max=10000, iter.max=10000, trace=ifelse(quiet==TRUE,0,1) ),
  #                        newtonsteps = newtonsteps,
  #                       ... )

  # Optimize
  results$opt = list( "par"=obj$par )
  for( i in seq_len(max(0,control$nlminb_loops)) ){
    if( isFALSE(control$quiet) ) message("Running nlminb_loop #", i)
    results$opt = nlminb( start = results$opt$par,
                  objective = obj$fn,
                  gradient = obj$gr,
                  control = list( eval.max = control$eval.max,
                                  iter.max = control$iter.max,
                                  trace = control$trace ) )
  }

  # Newtonsteps
  for( i in seq_len(max(0,control$newton_loops)) ){
    if( isFALSE(control$quiet) ) message("Running newton_loop #", i)
    g = as.numeric( obj$gr(results$opt$par) )
    h = optimHess(results$opt$par, fn=obj$fn, gr=obj$gr)
    results$opt$par = results$opt$par - solve(h, g)
    results$opt$objective = obj$fn(results$opt$par)
  }

  # Run sdreport
  if( isTRUE(control$getsd) ){
    if( isFALSE(control$quiet) ) message("Running sdreport")
    Hess_fixed = optimHess( par=results$opt$par, fn=obj$fn, gr=obj$gr )
    results$sdrep = sdreport( obj,
                              hessian.fixed=Hess_fixed,
                              getJointPrecision = control$getJointPrecision )
  }else{
    results$sdrep = NULL
  }

  results$report = obj$report()
  results$parhat = obj$env$parList()
  class(results) = "phylosem"
  return( results )
}

#' @title Detailed control for phylosem structure
#'
#' @description Define a list of control parameters.  Note that
#' the format of this input is likely to change more rapidly than that of
#' \code{\link{phylosem}}
#'
#' @param nlminb_loops Integer number of times to call \code{\link[stats]{nlminb}}.
#' @param newton_loops Integer number of Newton steps to do after running
#'        \code{\link[stats]{nlminb}}.
#' @param trace Parameter values are printed every `trace` iteration
#'        for the outer optimizer. Passed to
#'        `control` in \code{\link[stats]{nlminb}}.
#' @param eval.max Maximum number of evaluations of the objective function
#'        allowed. Passed to `control` in \code{\link[stats]{nlminb}}.
#' @param iter.max Maximum number of iterations allowed. Passed to `control` in
#'        \code{\link[stats]{nlminb}}.
#' @param getsd Boolean indicating whether to call \code{\link[TMB]{sdreport}}
#' @param run_model Boolean indicating whether to estimate parameters (the default), or
#'        instead to return the model inputs and compiled TMB object without running;
#' @param quiet Boolean indicating whether to run model printing messages to terminal or not;
#' @param getJointPrecision whether to get the joint precision matrix.  Passed
#'        to \code{\link[TMB]{sdreport}}.
#'
#' @return
#' An S3 object of class "phylosem_control" that specifies detailed model settings,
#' allowing user specification while also specifying default values
#'
#' @export
phylosem_control <-
function( nlminb_loops = 1,
          newton_loops = 1,
          trace = 0,
          eval.max = 1000,
          iter.max = 1000,
          getsd = TRUE,
          quiet = FALSE,
          run_model = TRUE,
          getJointPrecision = FALSE ){

  # Return
  structure( list(
    nlminb_loops = nlminb_loops,
    newton_loops = newton_loops,
    trace = trace,
    eval.max = eval.max,
    iter.max = iter.max,
    getsd = getsd,
    quiet = quiet,
    run_model = run_model,
    getJointPrecision = getJointPrecision
  ), class = "phylosem_control" )
}

#' @title Print parameter estimates and standard errors.
#'
#' @description Print parameter estimates
#' @param x Output from \code{\link{phylosem}}
#' @param ... Not used
#' @return prints (and invisibly returns) output from \code{\link[stats]{nlminb}}
#' @method print phylosem
#' @export
print.phylosem <- function(x, ...)
{
  cat("phylosem(.) result\n")
  if( "opt" %in% names(x) ){
    print( x$opt )
  }else{
    cat("`opt` not available in `phylosem`\n")
  }
  invisible(x$opt)
}

#' Extract path coefficients.
#'
#' @title Extract path coefficients
#' @param object Output from \code{\link{phylosem}}
#' @param standardized Whether to standardize regression coefficients
#' @param ... Not used
#' @return Data-frame listing all path coefficients, their parameter index and estimated values
#' @method coef phylosem
#' @export
coef.phylosem = function( object, standardized=FALSE, ... ){
  beta_z = object$opt$par[names(object$opt$par)=="beta_z"]
  RAM = object$obj$env$data$RAM
  if(nrow(RAM) != nrow(object$SEM_model)) stop("Check assumptions")
  for( i in which(RAM[,1]==1) ){
    if( standardized==TRUE ){
      beta_z[i] = beta_z[i] * abs(beta_z[which( RAM[,'from']==RAM[i,'from'] & RAM[,'to']==RAM[i,'from'] )])
      beta_z[i] = beta_z[i] / abs(beta_z[which( RAM[,'from']==RAM[i,'to'] & RAM[,'to']==RAM[i,'to'] )])
    }
  }
  # Report variances
  for( i in which(RAM[,1]==2) ){
    if( RAM[i,'from'] == RAM[i,'to'] ){
      beta_z[i] = beta_z[i]^2
    }
  }
  SEM_params = beta_z[ifelse(RAM[,4]==0, NA, RAM[,4])]
  SEM_params = ifelse( is.na(SEM_params), as.numeric(object$SEM_model[,3]), SEM_params )
  return( data.frame(Path=object$SEM_model[,1], Parameter=object$SEM_model[,2], Estimate=SEM_params ) )
}

#' @title Extract Variance-Covariance Matrix
#'
#' @description extract the covariance of fixed effects, or both fixed and random effects.
#'
#' @param object output from \code{phylosem}
#' @param which whether to extract the covariance among fixed effects, random effects, or both
#' @param ... ignored, for method compatibility
#' @importFrom stats vcov
#' @method vcov phylosem
#' @export
vcov.phylosem <-
function( object,
          which = c("fixed", "random", "both"),
          ...) {

  which = match.arg(which)

  if( which=="fixed" ){
    V = object$sdrep$cov.fixed
    if(is.null(V)){
      warning("Please re-run `phylosem` with `getsd=TRUE`, or confirm that the model is converged")
    }
  }
  if( which=="random" ){
    V = solve(object$obj$env$spHess(random=TRUE))
  }
  if( which=="both" ){
    H = object$sdrep$jointPrecision
    if(is.null(H)){
      warning("Please re-run `phylosem` with `getsd=TRUE` and `getJointPrecision=TRUE`, or confirm that the model is converged")
      V = NULL
    }else{
      V = solve(H)
    }
  }

  return( V )
}

# @title Extract the (marginal) log-likelihood of a phylosem model
#
# @return object of class \code{logLik} with attributes
#   \item{val}{log-likelihood}
#   \item{df}{number of parameters}
#' @importFrom stats logLik
#' @export
logLik.phylosem <- function(object, ...) {
  val = -1 * object$opt$objective
  df = length( object$opt$par )
  out = structure( val,
             df = df,
             class = "logLik")
  return(out)
}

#' @title summarize phylosem
#'
#' @description Summarize phylosem output from phylosem, including calculating intercepts at the tree root
#'
#' @param object Output from \code{\link{phylosem}}
#' @param ... Not used
#' @return Data-frame containing all estimated intercepts, path coefficients, and variance-covariance parameters
#'         as well as their standard errors
#' @method summary phylosem
#' @export
summary.phylosem = function( object, ... ){

  # Errors
  if(is.null(object$sdrep)) stop("Please re-run with `getsd=TRUE`")

  # Easy of use
  RAM = object$obj$env$data$RAM

  # Intercepts
  Intercepts = data.frame(
    Path = NA,
    VarName = paste0("Intercept_", colnames(object$data) ),
    Estimate = as.list(object$sdrep, "Estimate", report=TRUE)$intercept_j,
    StdErr = as.list(object$sdrep, "Std. Error", report=TRUE)$intercept_j
  )
  #rownames(Intercepts) = paste0("Intercept_", colnames(object$data) )

  # Slopes
  Slopes = data.frame(
    Path = object$SEM_model[which(RAM[,1]==1),1],
    VarName = object$SEM_model[which(RAM[,1]==1),2],
    Estimate = c(NA,as.list(object$sdrep, "Estimate")$beta_z)[RAM[which(RAM[,1]==1),4]+1],
    StdErr = c(NA,as.list(object$sdrep, "Std. Error")$beta_z)[RAM[which(RAM[,1]==1),4]+1]
  )
  # Plug in if fixed
  #Slopes$Estimate = ifelse( is.na(object$SEM_model[which(RAM[,1]==1),3]), Slopes$Estimate, as.numeric(object$SEM_model[which(RAM[,1]==1),3]) )
  Slopes$Estimate = ifelse( is.na(Slopes$Estimate), as.numeric(object$SEM_model[which(RAM[,1]==1),3]), Slopes$Estimate )
  # Unknown junk
  #rownames = object$SEM_model[which(RAM[,1]==1),2]
  #rownames( Slopes ) = rownames # ifelse( is.na(rownames), "TURNED OFF", rownames )

  # Covariances
  Variances = data.frame(
    Path = object$SEM_model[which(RAM[,1]==2),1],
    VarName = object$SEM_model[which(RAM[,1]==2),2],
    Estimate = c(NA,as.list(object$sdrep, "Estimate")$beta_z)[RAM[which(RAM[,1]==2),4]+1],
    StdErr = c(NA,as.list(object$sdrep, "Std. Error")$beta_z)[RAM[which(RAM[,1]==2),4]+1]
  )
  # Plug in if fixed
  #Variances$Estimate = ifelse( is.na(object$SEM_model[which(RAM[,1]==2),3]), Variances$Estimate, as.numeric(object$SEM_model[which(RAM[,1]==2),3]) )
  Variances$Estimate = ifelse( is.na(Variances$Estimate), as.numeric(object$SEM_model[which(RAM[,1]==2),3]), Variances$Estimate )
  # Unknown junk
  #rownames = object$SEM_model[which(RAM[,1]==2),1]
  #rownames( Variances ) = ifelse( is.na(rownames), "TURNED OFF", rownames )

  #
  Coefs = rbind( Intercepts, Slopes, Variances )
  Coefs = cbind( Coefs, "t.value"=abs(Coefs$Estimate)/Coefs$StdErr )
  Coefs = cbind( Coefs, "p.value"=2*(1-pnorm(abs(Coefs$t.value))) )

  out = list(
    call = object$call,
    coefficients = Coefs
  )

  # Return stuff
  class(out) = "summary.phylosem"
  return(out)
}

#' @title predict values for new tip
#'
#' @description Predict phylosem
#'
#' @inheritParams phytools::bind.tip
#'
#' @param x Output from \code{\link{phylosem}}
#' @param ... passed to \code{\link[phytools]{bind.tip}}
#'
#' @return predict new values
#'
#' @method predict phylosem
#' @export
#predict.phylosem <-
#function( x,
#          tip.label = "new_tip",
#          edge.length = NULL,
#          where = NULL,
#          ... ){
#
#  if( is.character(where) ){
#    where = which( c(x$tree$tip.label,x$tree$node.label) == where )
#    if(length(where)!=1) stop("`where` not found in `tree$tip.label` or `tree$node.label`")
#  }
#
#  # Return existing prediction OR rebuild
#  if( where <= Ntip(x$tree) ){
#    out = x$report$x_vj[where,]
#  }else{
#    # Add default edge.length
#    if(is.null(edge.length)){
#      if(is.ultrametric(x$tree)){
#        # node.depth, node.height, node.depth.edgelength
#        node_depth = node.depth.edgelength(x$tree)
#        edge_length = x$tree$edge.length
#        node_edgeout_length = edge_length[ match(1:(Ntip(x$tree)+Nnode(x$tree)),x$tree$edge[,1]) ]
#        node_and_edgeout_depth = node_depth + ifelse( is.na(node_edgeout_length), 0, node_edgeout_length )
#        tree_depth = max(node_and_edgeout_depth)
#        edge.length = tree_depth - node_depth[where]
#      }else{
#        stop("Must supply `edge.length`")
#      }
#    }
#
#    # build new tree
#    tree_new = phytools::bind.tip( x$tree,
#                                   tip.label = tip.label,
#                                   edge.length = edge.length,
#                                   where = where,
#                                   position = 0,
#                                   ... )
#
#    # refit
#    Args = as.list(x$call[-1])
#    Args$tree = tree_new
#    Args$startpar = x$opt$par
#    Args$quiet = TRUE
#    psem_new = do.call("phylosem", Args )
#
#    # extract
#    out = psem_new$report$x_vj[ which(psem_new$tree$tip.label == tip.label), ]
#  }
#
#  names(out) = colnames(x$data)
#  return( out )
#}
