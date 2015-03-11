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
test_that ('calcDist([basic]) works', {
  tree1 <- compute.brlen (stree (10, 'left'))
  tree2 <- compute.brlen (stree (10, 'left'))
  distances <- calcDist (tree1, tree2)
  # trees are identical
  expect_that (0, equals (distances[['PH85']]))
  expect_that (0, equals (distances[['score']]))
  expect_that (0, equals (distances[['dmat']]))
  # rearrange to be maximally distant
  tree1$tip.label <- tree1$tip.label[c (10,2,8,4,6,5,7,3,9,1)]
  distances <- calcDist (tree1, tree2)
  expect_that (1, equals (distances[['PH85']]))
  expect_that (1, equals (distances[['score']]))
  expect_that (0.6516854, equals (distances[['dmat']]))  # max for 10 tips
})

test_that ('mapNames([basic]) works', {
  # names with some typos and species not in tree
  names <- c ('Gorila gorila', 'Allenopithecus nigrowiridiss', 'Macaca madeppy',
              'Pan paniscus', 'Hylobates pileatus', 'Pygathrix bieti', 'Pan troglodytes',
              'Hylobates madeuppy', 'Pongo pingu', 'Homo neanderthelensis', 'Homo sapins')
  res <- mapNames (catarrhines, names, stats=TRUE)
  expect_that ('phylo', equals (class (res$tree)))
  # we expect the names to be resolved to the genus level and above
  expect_more_than (res$stat$mean.rank, 29)
})

test_that ('.searchTreeNames([basic]) works', {
  
})

test_that ('.addTip([basic]) works', {
  
})

test_that ('.extract([basic]) works', {
  
})