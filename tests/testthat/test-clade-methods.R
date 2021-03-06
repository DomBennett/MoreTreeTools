# Test clade-methods
# D.J. Bennett
# 24/06/2015

# LIBS
library (MoreTreeTools)
library (testthat)

# TEST DATA
data ('testclades')

# TEST FUNCTIONS
genTestClades <- function (n.clades=100, n.time=100, genfunc) {
  # generate a clade and plot with clade plot
  # genfunc can be random in different ways to test
  # the effectiveness of the plotting function
  # e.g. random uniform, random arching, random arching at different scales
  initclade <- rep (0, n.time)
  test.clades <- data.frame (initclade)
  starts <- floor (runif (min=2, max=n.time, n=n.clades))
  ends <- ceiling (runif (min=starts, max=n.time, n=n.clades))
  for (i in 1:n.clades) {
    temp <- initclade
    temp[starts[i]:ends[i]] <-
      genfunc ((ends[i] - starts[i]) + 1)
    test.clades <- cbind (test.clades, temp)
  }
  test.clades <- test.clades[,-1]
  test.clades
}

genUnif <- function (n) {
  # generate random uniform dist
  res <- ceiling (abs (rnorm (n)))
  res ^ max (res)
}

genArch <- function (n) {
  # create an arch on different scales
  half.n <- ceiling (n /2)
  first <- ceiling (abs (rnorm (n=half.n)))
  second <- ceiling (abs (rnorm (n=half.n)))
  res <- c (first[order (first, decreasing=FALSE)],
            second[order (second, decreasing=TRUE)])
  if (length (res) > n) {
    res <- res[-1]
  }
  res ^ max (res)
}

# RUNNING TESTS
context ('Testing \'calc-methods\'')
test_that ('getCladeSuccess([basic]) works', {
  
})

test_that ('calcCladeStats([basic]) works', {
  # calc stats for test clades
  # make sure CM and CG are as expected
  arch.stats <- calcCladeStats (arch.clades)
  arch.cm <- signif (mean (arch.stats$cm), 1)
  arch.cg <- signif (mean (arch.stats$cg), 1)
  expect_that (arch.cm, equals (0.5))
  unif.stats <- calcCladeStats (unif.clades)
  unif.cm <- signif (mean (unif.stats$cm), 1)
  unif.cg <- signif (mean (unif.stats$cg), 1)
  expect_that (unif.cm, equals (0.5))
  expect_that (unif.cg, is_more_than (arch.cg))
})

test_that ('plotClades([merge=FALSE]) works', {
  # simple class test
  p <- plotClades (clades=arch.clades)
  class.bool <- all (class (p) %in% c ("gg", "ggplot"))
  expect_true (class.bool)
})

test_that ('plotClades([merge=TRUE]) works', {
  # simple class test
  test.clades <- genTestClades (n.clades=10, n.time=10, genUnif)
  p <- plotClades (clades=test.clades, cids=1:ncol(test.clades), merge=TRUE)
  class.bool <- all (class (p) %in% c ("gg", "ggplot"))
  expect_true (class.bool)
  test.clades <- genTestClades (n.clades=10, n.time=10, genArch)
  p <- plotClades (clades=test.clades, cid=1:ncol(test.clades), merge=TRUE)
  class.bool <- all (class (p) %in% c ("gg", "ggplot"))
  expect_true (class.bool)
})

# unif.clades <- genTestClades (n.clades=100, n.time=100, genUnif)
# arch.clades <- genTestClades (n.clades=100, n.time=100, genArch)
# save (unif.clades, arch.clades, file='data/testclades.rda')
