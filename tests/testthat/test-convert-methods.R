# Test convert-methods

# LIBS
library(MoreTreeTools)
library(testthat)

# RUNNING
context ('Testing \'convert-methods\'')
test_that ('setAs(from=TreeMan, to=phylo) works', {
  tree <- treeman::randTree(10)
  tree <- as(tree, 'phylo')
  expect_true(is(tree, 'phylo'))
})

test_that ('setAs(from=phylo, to=TreeMan) works', {
  tree <- ape::rtree(10)
  tree <- as(tree, 'TreeMan')
  expect_true(is(tree, 'TreeMan'))
})