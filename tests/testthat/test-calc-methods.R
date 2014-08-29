## Test calc-methods
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

test_that ('calcED(type = \'FP\') works', {
  res <- calcED (catarrhines)
  # orangutan is the most distinct
  most.distinct <-
    as.character (rownames (res)[which (res[ ,1] == max (res[ ,1]))])
  expect_that ("Pongo pygmaeus", equals (most.distinct))
})

test_that ('calcED(type = \'PE\') works', {
  res <- calcED (catarrhines, type = 'PE')
  # orangutan is still the most distinct
  most.distinct <-
    as.character (rownames (res)[which (res[ ,1] == max (res[ ,1]))])
  expect_that ("Pongo pygmaeus", equals (most.distinct))
})

test_that ('calcED(type = \'ES\') works', {
  res <- calcED (catarrhines, type = 'ES')
  # orangutan is still the most distinct
  most.distinct <-
    as.character (rownames (res)[which (res[ ,1] == max (res[ ,1]))])
  expect_that ("Pongo pygmaeus", equals (most.distinct))
})