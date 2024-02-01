## ----include = FALSE----------------------------------------------------------
have_packages = all(sapply( c("phylopath"), FUN=requireNamespace))
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = have_packages  # https://r-pkgs.org/vignettes.html#sec-vignettes-eval-option
)
# To isntall
#   devtools::install_local( R'(C:\Users\James.Thorson\Desktop\Git\phylosem)', force=TRUE )
# Test build:
#  setwd( R'(C:\Users\James.Thorson\Desktop\Git\phylosem\vignettes)' ); devtools::build_rmd( R'(vignettes\phylopath.Rmd)' )
# Render example
#  library(rmarkdown); setwd( R'(C:\Users\James.Thorson\Desktop\Git\phylosem\vignettes)' ); render( "phylopath.Rmd", pdf_document())

## ----package_warning, include=isFALSE(have_packages)--------------------------
print("Must install ggplot2, phylopath, phylolm, ape")

## ----setup--------------------------------------------------------------------
data(rhino, rhino_tree, package = 'phylopath')

## -----------------------------------------------------------------------------
rhino_std <- rhino[-1]
rhino_std[] <- lapply(rhino_std, scale)

## -----------------------------------------------------------------------------
models_pp <- phylopath::define_model_set(
  one   = c(RS ~ DD),
  two   = c(DD ~ NL, RS ~ LS + DD),
  three = c(RS ~ NL),
  four  = c(RS ~ BM + NL),
  five  = c(RS ~ BM + NL + DD),
  six   = c(NL ~ RS, RS ~ BM),
  seven = c(NL ~ RS, RS ~ LS + BM),
  eight = c(NL ~ RS),
  nine  = c(NL ~ RS, RS ~ LS),
  .common = c(LS ~ BM, NL ~ BM, DD ~ NL)
)

## -----------------------------------------------------------------------------
models_ps <- list(
  one   = 'RS = b1 * DD',
  two   = 'DD = b1 * NL; RS = b2 * LS + b3 * DD',
  three = 'RS = b1 * NL',
  four  = 'RS = b * BM + b2 * NL',
  five  = 'RS = b1 * BM + b2 * NL + b3 * DD',
  six   = 'NL = b1 * RS; RS = b2 * BM',
  seven = 'NL = b1 * RS; RS = b2 * LS + b3 * BM',
  eight = 'NL = b1 * RS',
  nine  = 'NL = b1 * RS; RS = b2 * LS'
)
# we add the .common paths, by pasting them at the end of each of the model strings, e.g.:
models_ps <- lapply(
  models_ps, 
  \(x) paste(x, c('LS = b1_ * BM; NL = b2_ * BM; DD = b3_ * NL'), sep = '; ')
)

## ----eval=have_packages, results=FALSE, message=FALSE-------------------------
result_pp <- phylopath::phylo_path(
  models_pp, data = rhino_std, tree = rhino_tree, model = 'BM'
)

library(phylosem)
result_ps <- phylosem::compare_phylosem(
  models_ps, tree = rhino_tree, data = rhino_std
)

## -----------------------------------------------------------------------------
summary(result_pp)

## -----------------------------------------------------------------------------
sapply(result_ps, AIC) |> sort()

## -----------------------------------------------------------------------------
best_pp <- phylopath::best(result_pp)
best_ps <- phylosem::best(result_ps)

## ----eval=have_packages, fig.dim = c(6, 6), out.width = "50%"-----------------
plot(best_pp)
plot(as_fitted_DAG(best_ps))

## ----eval=have_packages, results=FALSE, message=FALSE, fig.dim = c(6, 6), out.width = "50%"----
phylopath::est_DAG(
  models_pp$five, data = rhino_std, tree = rhino_tree, model = 'lambda'
) |> plot()
phylosem::phylosem(
  models_ps$five, tree = rhino_tree, data = rhino_std, estimate_lambda = TRUE
) |> as_fitted_DAG() |> plot()

