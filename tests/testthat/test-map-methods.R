# Test map-methods
# D.J. Bennett
# 19/03/2015

# LIBS
library (MoreTreeTools)
library (testthat)

# TEST DATA
data ('catarrhines')

# RUNNING
context ('Testing \'map-methods\'')
test_that ('mapNames([basic]) works', {
  # names with some typos and species not in tree
  names <- c ('Gorila gorila', 'Allenopithecus nigrowiridiss', 'Macaca madeppy',
              'Pan paniscus', 'Hylobates pileatus', 'Pygathrix bieti', 'Pan troglodytes',
              'Hylobates madeuppy', 'Pongo pingu', 'Homo neanderthelensis', 'Homo sapins')
  res <- mapNames (catarrhines, names)
  # make sure it's a tree
  expect_that ('phylo', equals (class (res)))
  # make sure Homo n. has been placed on the pendant edge of Homo sapiens
  node <- getParent (catarrhines, tips = c ('Homo sapiens'))
  max.age <- getAge (catarrhines, node=node)
  node <- getParent (res, tips = c ('Homo sapins'))
  obs.age <- getAge (res, node=node)
  expect_less_than (obs.age, max.age)
})
test_that ('.mnMap([basic]) works', {
  paraenv <- new.env (parent=emptyenv ())
  paraenv$start.tree <- catarrhines
  paraenv$grow.tree <- catarrhines
  paraenv$datasource <- 4
  paraenv$matching.names <- c ('Pan troglodytes')
  paraenv$names <- c ('Pan troglodytes', 'Homo spiens', 'Gorilla grilla')
  sbjctenv <- new.env (parent=emptyenv ())
  MoreTreeTools:::.mnResolveUpdate (paraenv=paraenv, sbjctenv=sbjctenv)
  qrylist <- MoreTreeTools:::.mnResolve (names=c ('Homo spiens', 'Gorilla grilla'),
                                         paraenv)
  sbjctenv$resolved <- rbind (sbjctenv$resolved, qrylist$resolved)
  sbjctenv$lineages <- c (sbjctenv$lineages, qrylist$lineages)
  resenv <- new.env (parent=emptyenv ())
  resenv$trees <- list ()
  MoreTreeTools:::.mnMap (resenv=resenv, qrylist=qrylist, sbjctenv=sbjctenv,
                          paraenv=paraenv)
  coph <- cophenetic.phylo (resenv$trees[[1]])
  # expect chimp to be closer to human than gorilla
  expect_more_than (coph['Gorilla grilla', 'Homo spiens'],
                    coph['Pan troglodytes', 'Homo spiens'])
})
test_that ('.mnResolve([basic]) works', {
  paraenv <- list ()
  paraenv$datasource <- 4
  res <- MoreTreeTools:::.mnResolve (names=c ('Homo sapiens'), paraenv)
  expect_that (class (res$resolved), equals ('data.frame'))
  expect_that (class (res$lineages), equals ('list'))
})
test_that ('.mnResolveUpdate([basic]) works', {
  # set up some environments
  paraenv <- new.env (parent=emptyenv())
  sbjctenv <- new.env (parent=emptyenv())
  paraenv$grow.tree <- catarrhines
  paraenv$matching.names <- sample (catarrhines$tip.label, 10)
  # run twice, make sure second has more results
  MoreTreeTools:::.mnResolveUpdate (paraenv, sbjctenv)
  res1 <- nrow (sbjctenv$resolved)
  MoreTreeTools:::.mnResolveUpdate (paraenv, sbjctenv)
  res2 <- nrow (sbjctenv$resolved)
  expect_more_than (res2, res1)
})
test_that ('.mnSample([basic]) works', {
  # create a paraenv with tree
  paraenv <- list ()
  paraenv$grow.tree <- compute.brlen (rtree (10))
  paraenv$matching.names <- paraenv$grow.tree$tip.label
  # choose node up from root
  node <- getSize (paraenv$grow.tree) + 2
  # get all of its children and make already seen
  paraenv$deja.vues <- getChildren (paraenv$grow.tree, node)
  # names shouldn't have any deja.vues
  names <- MoreTreeTools:::.mnSample (paraenv)
  expect_false (any (names %in% paraenv$deja.vues))
})
test_that ('.mnTemporise([basic]) works', {
  paraenv <- list ()
  paraenv$datasource <- 4
  record <- MoreTreeTools:::.mnResolve (names=c ('Homo sapiens', 'Gallus gallus'),
                                        paraenv)
  res <- MoreTreeTools:::.mnTemporise (record=record, tree=catarrhines)
  # record was originally 2 rows
  expect_that (nrow (record$resolved), equals (2))
  # but temporise removes names not in tree
  expect_that (nrow (res$resolved), equals (1))
})
test_that ('.mnAddTip([basic]) works', {
  tree <- compute.brlen (rtree (10))
  tree <- MoreTreeTools:::.mnAddTip (tree, tip.i=10, new.name='t11')
  expect_that (tree$tip.label[11], equals ('t11'))
})
test_that ('.mnEarlyReturn([basic]) works', {
  sample.names <- sample (catarrhines$tip.label, 10)
  res <- MoreTreeTools:::.mnEarlyReturn (tree=catarrhines,
                                         names=sample.names,
                                         iterations=2)
  expect_that ('multiPhylo', equals (class (res)))
  expect_that (getSize (res[[1]]), equals (10))
})
test_that ('.mnExtract([basic]) works', {
  sample.names <- sample (catarrhines$tip.label, 10)
  res <- MoreTreeTools:::.mnExtract (tree=catarrhines,
                                     names=sample.names)
  expect_that ('phylo', equals (class (res)))
  expect_that (getSize (res), equals (10))
})