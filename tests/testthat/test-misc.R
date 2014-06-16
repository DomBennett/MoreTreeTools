## Test api tools
## D.J. Bennett
## 13/06/2014

## Libraries
library (MoreTreeTools)
library (testthat)

## Test data
test.lineages <- 
  list (c ('kingdomA', 'phylumA', 'classA',
           'orderA', 'familyA', 'genusA',
           'speciesA'),
        c ('kingdomA', 'phylumA', 'classA',
           'orderA', 'familyA', 'genusA',
           'speciesB'),
        c ('kingdomA', 'phylumA', 'classA',
           'orderA', 'familyA', 'genusA',
           'speciesC'),
        c ('kingdomA', 'phylumA', 'classA',
           'orderA', 'familyA', 'genusA',
           'speciesD'),
        c ('kingdomA', 'phylumA', 'classA',
           'orderA', 'familyA', 'genusB',
           'speciesE'),
        c ('kingdomA', 'phylumA', 'classA',
           'orderA', 'familyA', 'genusB',
           'speciesF'),
        c ('kingdomA', 'phylumA', 'classA',
           'orderA', 'familyA', 'genusC',
           'speciesG'),
        c ('kingdomA', 'phylumA', 'classA',
           'orderA', 'familyB', 'genusD',
           'speciesH'),
        c ('kingdomA', 'phylumA', 'classA',
           'orderA', 'familyC', 'genusE',
           'speciesI'),
        c ('kingdomA', 'phylumA', 'classA',
           'orderB', 'familyD', 'genusF',
           'speciesJ'))

## Running tests
context ('Testing \'misc\'')
test_that ('.findClade([basic]) works', {
  # class A is shared by all test species
  expect_that (MoreTreeTools:::.findClade (
    test.lineages), equals ('classA'))
})