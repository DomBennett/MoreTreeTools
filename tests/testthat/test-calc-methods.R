## Test api tools
## D.J. Bennett
## 13/06/2014

## Libraries
library (MoreTreeTools)
library (testthat)

## Test data
data ('catarrhines')

## Running tests
context ('Testing \'calc-methods\'')
test_that ('reduceTree([basic]) works', {
  # There are 23 Catarrhine genera
  res <- reduceTree (catarrhines, 'genus')
  expect_that (length (res$tip.label), equals (23))
})