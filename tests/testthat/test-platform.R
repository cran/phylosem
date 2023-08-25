

context("Testing cross platform and R version compatibility")

# Eastern Bering Sea pollcok
test_that("phylosem example is working ", {
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
          tree = rhino_tree )
  # Check objective function
  expect_equal( as.numeric(psem$opt$obj), 1087.686, tolerance=1e-2 )

  # Convert and plot using phylopath
  as_fitted_DAG(psem)
  as_sem(psem)
  as_phylo4d(psem)
})

