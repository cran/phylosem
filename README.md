## Phylogenetic structural equation models
Package _phylosem_ combines features from structural equation models (SEM), phylogenetic comparative methods (PCM), and generalized linear mixed models (GLMM).  By doing so, it incorporates a broad feature-set:

* Comparing multiple evolutionary models similar to _phylopath_
* Estimating trade-offs among multiple traits similar to _phylolm_
* Predicting missing trait values (and associated standard errors) similar to _Rphylopars_
* Estimating tradeoffs including recursive (cyclic) dependencies similar to package _sem_
* Applying ordination to multiple traits, similar to phylogenetic factor analysis in package _FishLife_

_phylosem_ is specifically intended as a minimal implementation, and uses standard packages for input/output formatting:

* Input: phylogenetic relatedness defined using class _phylo_ in package _ape_
* Input: structural trade-offs specified using syntax defined by package _sem_
* Output: visualizing trade-offs using _semPlot_, _diagrammeR_, and _ggraph_
* Output: assembling trait predictions and standard errors using _phylobase_
* Output: plotting trait predictions using _phylosignal_

Please see package vignettes for more details regarding syntax and features.

[![](https://cranlogs.r-pkg.org/badges/phylosem)](https://cran.r-project.org/package=phylosem)
[![](https://cranlogs.r-pkg.org/badges/grand-total/phylosem)](https://cran.r-project.org/package=phylosem)

[![DOI](https://zenodo.org/badge/534386257.svg)](https://zenodo.org/badge/latestdoi/534386257)



