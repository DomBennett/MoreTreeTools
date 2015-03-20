# Test server-methods
# D.J. Bennett
# 11/06/2014

# LIBS
library (MoreTreeTools)
library (testthat)

# TEST DATA
data ('hominoids')
hominoids <- hominoids$tip.label

# RUNNING
context ('Testing \'server-tools\'')
test_that ('.safeFromJSON([basic]) works', {
  # simply show that an error is thrown and handled
  expect_that (
    MoreTreeTools:::.safeFromJSON (
      url = 'dummyaddress', max.trys = 0), throws_error ())
})

test_that ('taxaResolve([basic]) works', {
  test.names <- c (hominoids, 'thisisnotaname')
  expected.dimensions <- c (17, 10)
  res <- MoreTreeTools:::taxaResolve (test.names)
  res <- res[complete.cases (res), ]
  expect_that (dim (res), equals (expected.dimensions))
})