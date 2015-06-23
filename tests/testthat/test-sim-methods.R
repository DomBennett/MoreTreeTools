# Test sim-methods
# D.J. Bennett
# 24/06/2015

# LIBS
library (MoreTreeTools)
library (testthat)

# RUNNING
context ('Testing \'sim-methods\'')
test_that ('.seedTree([basic]) works ...', {
  tree <- .seedTree (10, 10)
  # should have ten tips
  expect_that (length (tree$tip.label), equals (10))
  # should be 10 time units old
  expect_that (getAge (tree, node=length (tree$tip.label) + 1)$age,
               equals (10))
})
test_that ('runEDBMM([basic]) works ...', {
  # grow a tree by 10 tips
  tree <- runEDBMM (birth=2, death=1, stop.at=10, fossils=TRUE)
  # it should have 12 tips (because a default seed of two starts it)
  expect_that (getSize (drop.extinct(tree)), equals (12))
  # it should have extinct tips
  cat (paste0 ('NOTE! There is a 0.017 probability of \'runEDBMM([basic])\' test failing.
               Try running again if it does fail.'))
  expect_that (is.ultrametric (tree), is_false ())
})
test_that ('runEDBMM(record=TRUE) works ...', {
  # grow a tree by 20 tips, record tree growth every 10 tips added
  tree <- runEDBMM (birth=1, death=0.5, stop.at=20, sample.at=10,
                    record=TRUE)
  expect_that (getSize (tree[[1]]), equals (2))
  expect_that (getSize (tree[[2]]), equals (12))
  expect_that (getSize (tree[[3]]), equals (22))
})