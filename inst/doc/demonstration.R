## ----include = FALSE----------------------------------------------------------
have_packages = all(sapply( c("ggplot2","phylopath","phylolm","ape"), FUN=requireNamespace))
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = have_packages  # https://r-pkgs.org/vignettes.html#sec-vignettes-eval-option
)
# Test build:
#  devtools::build_rmd("vignettes/comparison.Rmd")
# Render example
#  library(rmarkdown); setwd( R'(C:\Users\James.Thorson\Desktop\Git\phylosem\vignettes)' ); render( "demonstration.Rmd", pdf_document())

## ----setup, echo=TRUE, warning=FALSE, message=FALSE---------------------------
library(phylosem)

## ----package_warning, echo=TRUE, eval=!have_packages--------------------------
#  message("Must install ggplot2, phylopath, phylolm, ape")

## ----eval=have_packages, echo=TRUE, message=FALSE, fig.width=6, fig.height=6----

# Compare using Pagel's kappa
library(phylopath)

# Run phylosem
model = "
  DD -> RS, p1
  BM -> LS, p2
  BM -> NL, p3
  NL -> DD, p4
"
psem = phylosem( sem = model,
          data = rhino[,c("BM","NL","DD","RS","LS")],
          estimate_ou = TRUE,
          estimate_lambda = TRUE,
          estimate_kappa = TRUE,
          tree = rhino_tree,
          getJointPrecision = TRUE,
          quiet = TRUE )

#
V = psem$opt$SD$cov.fixed
Rsub = cov2cor(V)[c('lnalpha','logitlambda','lnkappa'),c('lnalpha','logitlambda','lnkappa')]

knitr::kable(c("minimum_eigenvalue"=min(eigen(psem$opt$SD$jointPrecision)$values),
             "maximum_eigenvalue"=max(eigen(psem$opt$SD$jointPrecision)$values)), digits=3)
knitr::kable(Rsub, digits=3)

## ----eval=have_packages, echo=TRUE, message=FALSE, fig.width=6, fig.height=4----
library(ggplot2)
# Compile estimates and SEs
pdat = data.frame( "Estimate" = psem$opt$SD$par.fixed[c('lnalpha','logitlambda','lnkappa')],
              "StdErr" = sqrt(diag(V)[c('lnalpha','logitlambda','lnkappa')]) )
pdat = cbind( pdat, "Param" = rownames(pdat))

#
pdat$lower = pdat$Estimate - 1.96*pdat$StdErr
pdat$upper = pdat$Estimate + 1.96*pdat$StdErr
pdat$type = "estimated"

# Transform from log / logit-space to natural space
pdat2 = pdat
pdat2$Param = c("alpha", "lambda", "kappa")
pdat2['lnalpha',c("Estimate","lower","upper")] = exp(pdat2['lnalpha',c("Estimate","lower","upper")])
pdat2['lnkappa',c("Estimate","lower","upper")] = exp(pdat2['lnkappa',c("Estimate","lower","upper")])
pdat2['logitlambda',c("Estimate","lower","upper")] = plogis(as.numeric(pdat2['logitlambda',c("Estimate","lower","upper")]))
pdat2$type = "derived"

# Plot
ggplot( rbind(pdat,pdat2), aes( x=Param, y = Estimate, ymin = lower, ymax = upper)) +
  geom_pointrange(position = position_dodge(width = 0.6)) +
  theme_classic() +
  facet_grid( rows=vars(type), scales="free" )

## ----eval=have_packages, echo=TRUE, message=FALSE, fig.width=6, fig.height=6, out.width = "75%"----
# Settings
Ntree_config = c( 1e1, 1e2, 1e3, 1e4 )
Nreplicates = 5
sd_x = 0.3
sd_y = 0.3
b0_x = 1
b0_y = 0
b_xy = 1

# Simulate tree
set.seed(1)
Time_rcz = array(NA, dim=c(Nreplicates,length(Ntree_config),2), dimnames=list(NULL,"tree_size"=Ntree_config,"package"=c("phylolm","phylosem")) )

for( rI in seq_len(Nreplicates) ){
for( cI in seq_along(Ntree_config) ){

  # Simulate data
  tree = ape::rtree(n=Ntree_config[cI])
  x = b0_x + sd_x * phylolm::rTrait(n = 1, phy=tree)
  ybar = b0_y + b_xy*x
  y_normal = ybar + sd_y * phylolm::rTrait(n = 1, phy=tree)
  Data = data.frame(x=x, y=y_normal)[]

  # Run phylolm
  start_time = Sys.time()
  plm_bm = phylolm::phylolm(y ~ 1 + x, data=Data, phy=tree, model="BM" )
  Time_rcz[rI,cI,"phylolm"] = Sys.time() - start_time

  # Run phylosem
  start_time = Sys.time()
  psem_bm = phylosem( sem = "x -> y, p",
            data = Data,
            tree = tree,
            quiet = TRUE,
            newtonsteps = 0,
            getsd = FALSE )
  Time_rcz[rI,cI,"phylosem"] = Sys.time() - start_time
}}

# Format
df = apply( Time_rcz, MARGIN=2:3, FUN=mean )
df = cbind( expand.grid(dimnames(df)), "time_seconds"=as.vector(df) )

# Plot
library(ggplot2)
ggplot(data=df, aes(x=tree_size, y=time_seconds, group=package, color=package)) +
  geom_line() +
  geom_point() + scale_y_log10()

