#' @name runEDBMM
#' @title Evolutionary Distinctness Biased Markov Model (EDBMM)
#' @description Run simulation of tree growth using the Evolutionary
#' Distinctness Biased Markov Model (EDBMM).
#' @details None yet!
#' @template base_template
#' @param birth how many births per unit of branch length
#' @param death how many deaths per unit of branch length
#' @param stop.at what value to stop at
#' @param stop.by what to stop at time or taxa
#' @param seed.tree tree to start model with, else a random tree of 2 tips
#' @param bias what type of ED bias to use FP or PE
#' @param sig power modifier of ED speciation bias
#' @param eps power modifier of ED extinction bias
#' @param record return list of trees at specified sample points
#' @param fossils keep fossils in tree (T/F)
#' @param max.iterations maximum number of iterations to prevent inf looping
#' @param progress.bar plyr progress bar, record only
#' @export
#' @examples
#' # example.var <- exampleFun (test.data)

runEDBMM <- function (birth, death, stop.at, stop.by=c ('n', 't'),
                      seed.tree=NULL, bias=c ('FP', 'ES', 'PE'),
                      sig=-1, eps=-1, record=FALSE,
                      sample.at=stop.at*.1, fossils=FALSE,
                      max.iteration=10000, progress.bar='none') {
  # Wrapper function to allow recording of trees through time
  bias <- match.arg (bias)
  stop.by <- match.arg (stop.by)
  if (is.null (seed.tree)) {
    # create an arbitrary seed tree of size 2
    tree <- .seedTree (n = 2, age = (birth + death)/2)
  } else {
    # in this case we're going to build on top of the seed
    tree <- seed.tree
  }
  # create globals
  extinct <- NULL # vector of all extinct species
  max.node <- tree$Nnode
  max.tip <- getSize (tree)
  if (!record) {
    results <- .runEDBMM (birth, death, stop.at, tree,
                          bias, sig, eps, extinct, max.node,
                          max.tip, stop.by, fossils, max.iteration)
    return (results[['tree']])
  }
  iterations <- stop.at/sample.at # number of iterations
  trees <- list (tree) # collect trees
  runmodel <- function (i) {
    results <- .runEDBMM (birth, death, sample.at, tree,
                          bias, sig, eps, extinct, max.node,
                          max.tip, stop.by, fossils, max.iteration)
    tree <<- results[['tree']]
    max.node <<- results[['max.node']]
    max.tip <<- results[['max.tip']]
    extinct <<- results[['extinct']]
    trees <<- c (trees, list (tree))
  }
  m_ply (.data = (i = 1:iterations), .fun = runmodel,
         .progress = progress.bar)
  class (trees) <- 'multiPhylo'
  return (trees)
}

.runEDBMM <- function (birth, death, stop.at, seed.tree,
                       bias, sig, eps, extinct, max.node,
                       max.tip, stop.by, fossils, max.iteration) {
  # Grow tree using a markov model with diversification rates
  #  based on evolutionary distinctness (EDBMM)
  # Internal functions
  randomTip <- function (add = TRUE) {
    # Return a random tip based on bias and whether
    #  it is being added or not
    .calcED <- function (tree) {
      if (!is.null (extinct)) {
        extant.tree <- drop.tip (tree, tip = extinct)
        eds <- calcED (extant.tree, type = bias)
      } else {
        eds <- calcED (tree, type = bias)
      }
      eds
    }
    probs <- .calcED (tree)
    if (add) {
      # apply sigma if adding
      probs[ ,1] <- probs[ ,1]^sig
    } else {
      # apply epsilon if subtracting
      probs[ ,1] <- probs[ ,1]^eps
    }
    # return a species name based on probs
    sample (rownames (probs), size = 1, prob = probs[ ,1])
  }
  add <- function () {
    # add a tip to the tree
    extant.tips <- which (!tree$tip.label %in% extinct)
    # time passed is 1/length (extant.tips) -- so that
    #  births and deaths are scaled to 1 unit of branch length
    time.passed <- branch.growth/length (extant.tips)
    # find tip edges for all extant tips
    tip.edges <- which (tree$edge[ ,2] %in% extant.tips)
    # extant edges grow by 1/length (extant.tips)
    tree$edge.length[tip.edges] <-
      tree$edge.length[tip.edges] + time.passed
    # find random species to speciate
    to.speciate <- randomTip ()
    # find its tip edge
    to.speciate <- which (tree$tip.label == to.speciate)
    to.speciate <- which (tree$edge[ ,2] %in% to.speciate)
    # new node must have unique labels
    new.node.label <- paste0 ('n', max.node + 1)
    new.tip.label <- paste0 ('t', max.tip + 1)
    # new node.age is always time.passed -- one time step ago.
    node.age <- time.passed
    # add new tip at random edge
    tree <<- addTip (tree = tree, edge = to.speciate,
                     tip.name = new.tip.label,
                     node.age = node.age,
                     node.label = new.node.label)
    # add 1 to globals to keep names unique!
    max.node <<- max.node + 1
    max.tip <<- max.tip + 1
    # add to n and time
    t <<- t + time.passed
    n <<- n + 1
  }
  drop <- function () {
    # choose a species to drop
    to.drop <- as.character(randomTip (add = FALSE))
    if (fossils) {
      # either add to the extinct vector if fossils are to be kept
      extinct <<- c (extinct, to.drop)
    } else {
      # or drop the species
      tree <<- drop.tip (tree, to.drop)
    }
    n <<- n - 1
  }
  run <- function () {
    # run function
    # randomly add or drop
    add.bool <- sample (c (TRUE, FALSE), size = 1,
                        prob = c (birth, death))
    if (add.bool) {
      add ()
    } else {
      n.extant <- getSize (tree) - length (extinct)
      if (n.extant > 2) {
        # only drop species if there are more than
        # two species in a tree, this might create edge
        # effects for certain parameters
        drop ()
      }
    }
  }
  runForT <- function (i) {
    # run and stop after t
    run ()
    if (t >= stop.at) {
      stop (exp.message)
    }
  }
  runForN <- function (i) {
    # run and stop after n
    run ()
    if (n >= stop.at) {
      stop (exp.message)
    }
  }
  # set n and t to zero
  n <- t <- 0
  # calculate branch growth
  branch.growth <- 1/(birth + death)
  # set tree as seed.tree
  tree <- seed.tree
  # expected stop message
  exp.message <- 'Max reached'
  # run m_ply until stop ()
  if (stop.by == 't') {
    try (expr = m_ply (.data = (i = 1:max.iteration), .fun = runForT),
         silent = TRUE)
  } else {
    try (expr = m_ply (.data = (i = 1:max.iteration), .fun = runForN),
         silent = TRUE)
  }
  error.message <- geterrmessage ()
  # if last error is not exp.message raise it
  if (!grepl (exp.message, error.message)) {
    stop (error.message)
  }
  list (tree = tree, max.node = max.node, max.tip = max.tip,
        extinct = extinct)
}
.seedTree <- function (n, age) {
  # Create a seed tree for growing with n tips
  tree <- rtree (n)
  # add edge lengths
  tree <- compute.brlen (tree)
  tree.age <- getAge (tree, node = length (tree$tip.label) + 1)$age
  # re-calculate edge lengths based on requested tree age
  tree$edge.length <- tree$edge.length/(tree.age/age)
  # add labels
  tree$node.label <- paste0 ('n', 1:tree$Nnode)
  tree
}