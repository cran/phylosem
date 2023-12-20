## ----include = FALSE----------------------------------------------------------
have_packages = all(sapply( c("ggplot2","phylolm","phylopath","Rphylopars","phyr","TreeTools"), FUN=requireNamespace))
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = have_packages  # https://r-pkgs.org/vignettes.html#sec-vignettes-eval-option
)
# devtools::build_rmd("vignettes/comparison.Rmd") # to test build
# setwd( R'(C:\Users\James.Thorson\Desktop\Git\phylosem\vignettes)' ); rmarkdown::render( "comparison.Rmd", rmarkdown::pdf_document()) # to build PDF, and perhaps tinytex::latexmk() to install dependencies

## ----setup, echo=TRUE, warning=FALSE, message=FALSE---------------------------
library(phylosem)

## ----package_warning, include=!have_packages----------------------------------
message("Must install ggplot2, phylopath, phylolm, Rphylopars, TreeTools")

## ----eval=have_packages, echo=TRUE, message=FALSE, fig.width=6, fig.height=6----
# Settings
Ntree = 100
sd_x = 0.3
sd_y = 0.3
b0_x = 1
b0_y = 0
b_xy = 1

# Simulate tree
set.seed(1)
tree = ape::rtree(n=Ntree)

# Simulate data
x = b0_x + sd_x * phylolm::rTrait(n = 1, phy=tree)
ybar = b0_y + b_xy*x
y_normal = ybar + sd_y * phylolm::rTrait(n = 1, phy=tree)

# Construct, re-order, and reduce data
Data = data.frame(x=x,y=y_normal)[]

# Compare using BM model
start_time = Sys.time()
plm_bm = phylolm::phylolm(y ~ 1 + x, data=Data, phy=tree, model="BM" )
Sys.time() - start_time
knitr::kable(summary(plm_bm)$coefficients, digits=3)

start_time = Sys.time()
psem_bm = phylosem( sem = "x -> y, p",
          data = Data,
          tree = tree,
          quiet = TRUE )
Sys.time() - start_time
knitr::kable(summary(psem_bm)$coefficients, digits=3)

# Compare using OU
start_time = Sys.time()
plm_ou = phylolm::phylolm(y ~ 1 + x, data=Data, phy=tree, model="OUrandomRoot" )
Sys.time() - start_time

start_time = Sys.time()
psem_ou = phylosem( sem = "x -> y, p",
          data = Data,
          tree = tree,
          estimate_ou = TRUE,
          quiet = TRUE )
Sys.time() - start_time

knitr::kable(summary(psem_ou)$coefficients, digits=3)
knitr::kable(summary(plm_ou)$coefficients, digits=3)
knitr::kable(c( "phylolm_alpha"=plm_ou$optpar, 
                "phylosem_alpha"=exp(psem_ou$parhat$lnalpha) ), digits=3)

# Compare using Pagel's lambda
start_time = Sys.time()
plm_lambda = phylolm::phylolm(y ~ 1 + x, data=Data, phy=tree, model="lambda" )
Sys.time() - start_time

start_time = Sys.time()
psem_lambda = phylosem( sem = "x -> y, p",
          data = Data,
          tree = tree,
          estimate_lambda = TRUE,
          quiet = TRUE )
Sys.time() - start_time

knitr::kable(summary(psem_lambda)$coefficients, digits=3)
knitr::kable(summary(plm_lambda)$coefficients, digits=3)
knitr::kable(c( "phylolm_lambda"=plm_lambda$optpar, 
                "phylosem_lambda"=plogis(psem_lambda$parhat$logitlambda) ), digits=3)

# Compare using Pagel's kappa
start_time = Sys.time()
plm_kappa = phylolm::phylolm(y ~ 1 + x, data=Data, phy=tree, model="kappa", lower.bound = 0, upper.bound = 3 )
Sys.time() - start_time

start_time = Sys.time()
psem_kappa = phylosem( sem = "x -> y, p",
          data = Data,
          tree = tree,
          estimate_kappa = TRUE,
          quiet = TRUE )
Sys.time() - start_time

knitr::kable(summary(psem_kappa)$coefficients, digits=3)
knitr::kable(summary(plm_kappa)$coefficients, digits=3)
knitr::kable(c( "phylolm_kappa"=plm_kappa$optpar, 
                "phylosem_kappa"=exp(psem_kappa$parhat$lnkappa) ), digits=3)

## ----eval=have_packages, echo=TRUE, message=FALSE, fig.width=6, fig.height=6, out.width = "75%"----
# Settings
Ntree = 100
sd_x = 0.3
sd_y = 0.3
b0_x = 1
b0_y = 0
b_xy = 1

# Simulate tree
set.seed(1)
tree = ape::rtree(n=Ntree)

# Simulate data
x = b0_x + sd_x * phylolm::rTrait(n = 1, phy=tree)
ybar = b0_y + b_xy*x
y_normal = ybar + sd_y * phylolm::rTrait(n = 1, phy=tree)
y_pois = rpois( n=Ntree, lambda=exp(y_normal) )

# Construct, re-order, and reduce data
Data = data.frame(x=x,y=y_pois)

# Compare using phylolm::phyloglm
pglm = phylolm::phyloglm(y ~ 1 + x, data=Data, phy=tree, method="poisson_GEE" )
knitr::kable(summary(pglm)$coefficients, digits=3)

#
pglmm = phyr::pglmm_compare(
  y ~ 1 + x,
  family = "poisson",
  data = Data,
  phy = tree )
knitr::kable(summary(pglmm), digits=3)

#
pgsem = phylosem( sem = "x -> y, p",
          data = Data,
          family = c("fixed","poisson"),
          tree = tree,
          quiet = TRUE )
knitr::kable(summary(pgsem)$coefficients, digits=3)

## ----echo=TRUE, message=FALSE, fig.width=6, fig.height=6, out.width = "75%", eval=FALSE----
#  # Comare using Bayesian implementation in brms
#  library(brms)
#  Amat <- ape::vcv.phylo(tree)
#  Data$tips <- rownames(Data)
#  mcmc <- brm(
#    y ~ 1 + x + (1 | gr(tips, cov = A)),
#    data = Data, data2 = list(A = Amat),
#    family = 'poisson',
#    cores = 4
#  )
#  knitr::kable(fixef(mcmc), digits = 3)
#  
#  # Plot them together
#  library(ggplot2)
#  pdat <- rbind.data.frame(
#    coef(summary(pglm))[, 1:2],
#    data.frame(Estimate = pglmm$B, StdErr = pglmm$B.se),
#    setNames(as.data.frame(fixef(mcmc))[1:2], c('Estimate', 'StdErr')),
#    setNames(summary(pgsem)$coefficients[2:3, 3:4], c('Estimate', 'StdErr'))
#  )
#  pdat$Param <- rep(c('Intercept', 'Slope'), 4)
#  pdat$Method <- rep( c('phylolm::phyloglm', 'phyr::pglmm_compare',
#                       'brms::brm', 'phylosem::phylosem'), each = 2)
#  figure = ggplot(pdat, aes(
#    x = Estimate, xmin = Estimate - StdErr,
#    xmax = Estimate + StdErr, y = Param, color = Method
#  )) +
#    geom_pointrange(position = position_dodge(width = 0.6)) +
#    theme_classic() +
#    theme(panel.grid.major.x = element_line(), panel.grid.minor.x = element_line())

## ----include = FALSE, eval=FALSE----------------------------------------------
#  saveRDS( figure, file=file.path(R'(C:\Users\James.Thorson\Desktop\Git\phylosem\vignettes)',"brms.RDS") )
#  saveRDS( pdat, file=file.path(R'(C:\Users\James.Thorson\Desktop\Git\phylosem\vignettes)',"pdat.RDS") )

## ----eval=have_packages, message=FALSE, fig.width=6, fig.height=6, out.width = "75%", echo=FALSE----
library(ggplot2)
#figure = readRDS( file.path(system.file("brms",package="phylosem"),"brms.RDS") )
#figure
pdat = readRDS( file.path(system.file("brms",package="phylosem"),"pdat.RDS") )
ggplot(pdat, aes(
  x = Estimate, xmin = Estimate - StdErr,
  xmax = Estimate + StdErr, y = Param, color = Method
)) +
  geom_pointrange(position = position_dodge(width = 0.6)) +
  theme_classic() +
  theme(panel.grid.major.x = element_line(), panel.grid.minor.x = element_line())

## ----eval=have_packages, echo=TRUE, message=FALSE, fig.width=6, fig.height=6----
# Settings
Ntree = 100
sd_x = 0.3
sd_y = 0.3
b0_x = 1
b0_y = 0
b_xy = 1

# Simulate tree
set.seed(1)
tree = ape::rtree(n=Ntree)

# Simulate data
x = b0_x + sd_x * phylolm::rTrait(n = 1, phy=tree)
ybar = b0_y + b_xy*x
y_normal = ybar + sd_y * phylolm::rTrait(n = 1, phy=tree)
y_binom = rbinom( n=Ntree, size=1, prob=plogis(y_normal) )

# Construct, re-order, and reduce data
Data = data.frame(x=x,y=y_binom)

#
pglmm = phyr::pglmm_compare(
  y ~ 1 + x,
  family = "binomial",
  data = Data,
  phy = tree )
knitr::kable(summary(pglmm), digits=3)

#
pgsem = phylosem( sem = "x -> y, p",
          data = Data,
          family = c("fixed","binomial"),
          tree = tree,
          quiet = TRUE )
knitr::kable(summary(pgsem)$coefficients, digits=3)

## ----eval=have_packages, echo=TRUE, message=FALSE, fig.width=6, fig.height=6----
library(phylopath)
library(phylosem)

# make copy of data that's rescaled
rhino_scaled = rhino
  rhino_scaled[,c("BM","NL","LS","DD","RS")] = scale(rhino_scaled[,c("BM","NL","LS","DD","RS")])

# Fit and plot using phylopath
dag <- DAG(RS ~ DD, LS ~ BM, NL ~ BM, DD ~ NL)
start_time = Sys.time()
result <- est_DAG( DAG = dag,
                    data = rhino,
                    tree = rhino_tree,
                    model = "BM",
                    measurement_error = FALSE )
Sys.time() - start_time
plot(result)

# Fit and plot using phylosem
model = "
  DD -> RS, p1
  BM -> LS, p2
  BM -> NL, p3
  NL -> DD, p4
"
start_time = Sys.time()
psem = phylosem( sem = model,
          data = rhino_scaled[,c("BM","NL","DD","RS","LS")],
          tree = rhino_tree,
          quiet = TRUE )
Sys.time() - start_time
plot( as_fitted_DAG(psem) )

## ----eval=have_packages, echo=TRUE, message=FALSE, fig.width=6, fig.height=6----
library(sem)
library(TreeTools)

# Simulation parameters
n_obs = 50
# Intercepts
a1 = 1
a2 = 2
a3 = 3
a4 = 4
# Slopes
b12 = 0.3
b23 = 0
b34 = 0.3
# Standard deviations
s1 = 0.1
s2 = 0.2
s3 = 0.3
s4 = 0.4

# Simulate data
E1 = rnorm(n_obs, sd=s1)
E2 = rnorm(n_obs, sd=s2)
E3 = rnorm(n_obs, sd=s3)
E4 = rnorm(n_obs, sd=s4)
Y1 = a1 + E1
Y2 = a2 + b12*Y1 + E2
Y3 = a3 + b23*Y2 + E3
Y4 = a4 + b34*Y3 + E4
Data = data.frame(Y1=Y1, Y2=Y2, Y3=Y3, Y4=Y4)

# Specify path diagram (in this case, using correct structure)
equations = "
  Y2 = b12 * Y1
  Y4 = b34 * Y3
"
model <- specifyEquations(text=equations, exog.variances=TRUE, endog.variances=TRUE)

# Fit using package:sem
start_time = Sys.time()
Sem <- sem(model, data=Data)
Sys.time() - start_time

# Specify star phylogeny
tree_null = TreeTools::StarTree(n_obs)
  tree_null$edge.length = rep(1,nrow(tree_null$edge))
  rownames(Data) = tree_null$tip.label

# Fit using phylosem
start_time = Sys.time()
psem = phylosem( data = Data,
          sem = equations,
          tree = tree_null,
          quiet = TRUE )
Sys.time() - start_time

## ----eval=have_packages, echo=FALSE, results='asis'---------------------------
knitr::kable(coef(Sem, standardized=TRUE)[1:2], digits=3)
knitr::kable(coef(psem, standardized=TRUE)[1:2,], digits=3)

## ----eval=have_packages, echo=FALSE, results='asis'---------------------------
knitr::kable(coef(Sem, standardized=FALSE), digits=3)
knitr::kable(coef(psem, standardized=FALSE), digits=3)

## ----eval=have_packages, echo=TRUE, message=FALSE, fig.width=6, fig.height=6----
library(Rphylopars)

# Format data, within no values for species t1
Data = rhino[,c("BM","NL","DD","RS","LS")]
  rownames(Data) = tree$tip.label
Data['t1',] = NA

# fit using phylopars
start_time = Sys.time()
pars <- phylopars( trait_data = cbind(species=rownames(Data),Data),
                  tree = tree,
                  pheno_error = FALSE,
                  phylo_correlated = TRUE,
                  pheno_correlated = FALSE)
Sys.time() - start_time

# Display estimates for missing values
knitr::kable(cbind( "Estimate"=pars$anc_recon["t1",], "Var"=pars$anc_var["t1",] ), digits=3)

# fit using phylosem
start_time = Sys.time()
psem = phylosem( data = Data,
                 tree = tree,
                 sem = "",
                 covs = "BM, NL, DD, RS, LS",
                 quiet = TRUE )
Sys.time() - start_time

# Display estimates for missing values
knitr::kable(cbind(
  "Estimate"=as.list(psem$opt$SD,"Estimate")$x_vj[ match("t1",tree$tip.label), ],
  "Var"=as.list(psem$opt$SD,"Std. Error")$x_vj[ match("t1",tree$tip.label), ]^2
), digits=3)

