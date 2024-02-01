

context("Testing cross platform and R version compatibility")

# Eastern Bering Sea pollcok
test_that("phylosem example is working ", {
  #
  checkDepPackageVersion( warn=FALSE )

  #skip_on_ci()
  data(rhino, rhino_tree, package="phylopath")

  # Run phylosem
  model = "
    DD -> RS, p1
    BM -> LS, p2
    BM -> NL, p3
    NL -> DD, p4
  "
  psem = phylosem( sem = model,
          data = rhino[,c("BM","NL","DD","RS","LS")],
          tree = rhino_tree,
          control = phylosem_control( getJointPrecision=TRUE ) )
  # Check objective function
  expect_equal( as.numeric(psem$opt$obj), 1087.686, tolerance=1e-2 )

  # Convert and plot using phylopath
  as_fitted_DAG(psem)
  as_sem(psem)
  as_phylo4d(psem)

  #
  logLik(psem)
  AIC(psem)
  summary(psem)
  print(sem)
  coef(psem)
  vcov(psem, which="fixed")
  vcov(psem, which="random")
  vcov(psem, which="both")

  #
  psem_full = phylosem( sem = model,
        data = rhino[,c("BM","NL","DD","RS","LS")],
        tree = rhino_tree,
        estimate_ou = TRUE,
        estimate_lambda = TRUE,
        estimate_kappa = TRUE,
        control = phylosem_control( getsd = FALSE,
                                    quiet = TRUE ) )

  #
  model2 = "
    LS -> BM, p2
    BM -> NL, p3
    NL -> DD, p4
    DD -> RS, p1
  "
  myset = compare_phylosem( sem_set = list("one"=model,"two"=model2),
          data = rhino[,c("BM","NL","DD","RS","LS")],
          tree = rhino_tree )
  best(myset)
  choice(myset,1)
  average(myset)
})

