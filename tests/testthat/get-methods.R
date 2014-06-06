## Test get methods
## D.J. Bennett
## 05/05/2014

context ('getChildren')
test_that ('getChildren basic run',{
  data ('hominoids')
  great.apes <-
    c ("Pongo pygmaeus", "Gorilla gorilla",
       "Homo sapiens", "Pan paniscus",
       "Pan troglodytes")
  expect_that (getChildren (hominoids, 20),
               equals (great.apes))
  expect_that (getChildren (hominoids, '20'),
               throws_error ())
  expect_that (getChildren (hominoids, 31),
               throws_error ())
})