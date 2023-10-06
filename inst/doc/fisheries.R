## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
# To test build:
#   setwd( R'(C:\Users\James.Thorson\Desktop\Git\phylosem)' ); devtools::build_rmd("vignettes/fisheries.Rmd")
# To build PDF: 
#   library(rmarkdown); setwd( R'(C:\Users\James.Thorson\Desktop\Git\phylosem\vignettes)' ); render( "fisheries.Rmd", pdf_document()) 

## ----echo=TRUE, results='hide', message=FALSE, fig.width=6, fig.height=6------
# Load packages
library(phylosem)
library(fishtree)

# Download tree
out = fishtree_complete_phylogeny()
tree = out[[1]]

# Load data object
data( Mlifehist_ver1_0 )
Data = Mlifehist_ver1_0

# Reformat to match tree$tip.label
Data$Genus_species = factor( paste0(Data$Genus, "_", Data$Species) )

# Drop duplicates ... not dealing with variation among stocks within species
Data = Data[match(unique(Data$Genus_species),Data$Genus_species), ]

# log-transform to simplify later syuntax
Data = cbind( Data, "logM" = log(Data[,'M']),
                    "logK" = log(Data[,'K']),
                    "logtmax" = log(Data[,'tmax']),
                    "logLinf" = log(Data[,'Linf']) )

# Identify species in both datasets
species_to_use = intersect( tree$tip.label, Data$Genus_species )
species_to_drop = setdiff( Data$Genus_species, tree$tip.label )

# Drop tips not present in trait-data
# Not strictly necessary, but helpful to simplify later plots
tree = ape::keep.tip( tree, tip=species_to_use )

# Drop trait-data not in phylogeny
# Necessary to define correlation among data
rows_to_use = which( Data$Genus_species %in% species_to_use )
Data = Data[rows_to_use,]

# Only include modeled variables in trait-data passed to phylosem
rownames(Data) = Data$Genus_species
Data = Data[,c('logM','logK','logtmax','logLinf')]

## ----echo=TRUE, results='hide', message=FALSE, fig.width=6, fig.height=6------
# Specify SEM structure
sem_structure = "
  logK -> logtmax, b1
  logLinf -> logtmax, b2
  logtmax -> logM, a
"

# Grid-search model selection using AIC for transformations
Grid = expand.grid( "OU" = c(FALSE,TRUE),
                    "lambda" = c(FALSE,TRUE),
                    "kappa" = c(FALSE,TRUE) )
psem_grid = NULL
for( i in 1:nrow(Grid)){
  psem_grid[[i]] = phylosem( data=Data,
                   tree = tree,
                   sem = sem_structure,
                   estimate_ou = Grid[i,'OU'],
                   estimate_lambda = Grid[i,'lambda'],
                   estimate_kappa = Grid[i,'kappa'],
                   quiet = TRUE )
}

# Extract AIC for each model and rank-order by parsimony
Grid$AIC = sapply( psem_grid, \(m) m$opt$AIC )
Grid = Grid[order(Grid$AIC,decreasing=FALSE),]

# Select model with lowest AIC
psem_best = psem_grid[[as.numeric(rownames(Grid[1,]))]]

## ----echo=FALSE, message=FALSE, fig.width=6, fig.height=6---------------------
knitr::kable(Grid, digits=3, row.names=FALSE)
knitr::kable(summary(psem_best)$coefficients, digits=3)

## ----echo=TRUE, message=FALSE, fig.width=4, fig.height=4----------------------
# Plot path diagram
my_fitted_DAG = as_fitted_DAG(psem_best)
plot(my_fitted_DAG, type="color")

# Total, direct, and indirect effects
my_sem = as_sem(psem_best)
effects(my_sem)

## ----eval=FALSE, echo=TRUE, message=FALSE, fig.width=6, fig.height=8----------
#  # Load for plotting,
#  # requireNamespace("phylosignal",quietly=TRUE)
#  # https://r-pkgs.org/vignettes.html#sec-vignettes-eval-option
#  #library(phylosignal)
#  
#  # Plot using phylobase
#  #my_phylo4d = as_phylo4d( psem_best )
#  #barplot(my_phylo4d)

## ----echo=TRUE, results='hide', message=FALSE, fig.width=4, fig.height=4------
library(ape)
Data = Mlifehist_ver1_0

# Make taxonomic factors
Data$Genus_species = factor( paste0(Data$Genus, "_", Data$Species) )
Data$Genus = factor( Data$Genus )
Data$Family = factor( Data$Family )
Data$Order = factor( Data$Order )

# Make taxonomic tree
tree = ape::as.phylo( ~Order/Family/Genus/Genus_species, data=Data, collapse=FALSE)
tree$edge.length = rep(1,nrow(tree$edge))
tree = collapse.singles(tree)
tmp = root(tree, node=ape::Ntip(tree)+1 )

# Drop duplicates ... not dealing with variation among stocks within species
Data = Data[match(unique(Data$Genus_species),Data$Genus_species), ]

# log-transform to simplify later syuntax
Data = cbind( Data, "logM" = log(Data[,'M']),
                    "logK" = log(Data[,'K']),
                    "logtmax" = log(Data[,'tmax']),
                    "logLinf" = log(Data[,'Linf']) )

# Only include modeled variables in trait-data passed to phylosem
rownames(Data) = Data$Genus_species
Data = Data[,c('logM','logK','logtmax','logLinf')]

# Fit model
psem_taxon = phylosem( data=Data,
                 tree = tree,
                 sem = sem_structure,
                 estimate_ou = TRUE,
                 quiet = TRUE )

# Plot path diagram
my_fitted_DAG = as_fitted_DAG(psem_taxon)
plot(my_fitted_DAG, type="color")

