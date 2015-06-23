#' @name getCladeSuccess
#' @title Clade success through time
#' @description Return a dataframe of clade tip number from a list
#' of trees.
#' @details The function takes a list of trees and counts the success
#' of each clade, returning a dataframe of number of tips for each
#' tree for each named clade. It requires nodes to be labeled in
#' the 'node.label' slot of the Phylo class.
#' 
#' Takes list of trees generated from runEDBMM(record=TRUE).
#' 
#' Note, a clade in this function is defined as any node with two
#' descending branches each with extant tips.
#' @template base_template
#' @param trees
#' @param time.intervals vector of time represented by each tree in trees
#' @export
#' @examples
#' # example.var <- exampleFun (test.data)

getCladeSuccess <- function (trees, time.intervals=1:length(trees)) {
  # Return dataframe of clades through time from list of trees
  # with named nodes
  .get <- function (i) {
    .countChildren (trees[[i]])
  }
  clade.successes <- mlply (.data = data.frame (i = 1:length (trees)),
                            .fun = .get)
  .reformat (clade.successes, time.intervals)
}
.countChildren <- function (tree) {
  # Count the number of extant children for every node
  .count <- function (node.label) {
    node <- which (tree$node.label == node.label)
    node <- length (tree$tip.label) + node
    if (sum (tree$edge[ ,1] == node) > 1) {
      # if the node has two descending edges
      node.children <- getChildren (tree, node)
      extant <- node.children[node.children %in% extants]
      return (length (extant))
    }
    return (0)
  }
  # all internal nodes
  node.labels <- tree$node.label
  # make sure only extant clades are counted
  extants <- getExtant (tree)
  res <- mdply (.data = data.frame (node.label = node.labels),
                .fun = .count)
  colnames (res) <- c ('node', 'n.children')
  res
}
.reformat <- function (clade.performance, time.intervals) {
  # Take list of list of clade performances and convert
  #  to a dataframe.
  .getTime <- function (i, node) {
    # get success for node at a time point
    data <- clade.performance[[i]]
    if (any (data$node == node)) {
      return (data[data$node == node, 2])
    }
    0
  }
  .getNode <- function (node) {
    mdply (.data = data.frame (
      i = 1:length (clade.performance)),
      .fun = .getTime, node)[ ,2]
  }
  .addNode <- function (node) {
    # add success for node at all time points
    # for a res dataframe
    node <- as.character (node)
    node.success <- .getNode (node)
    res[node] <- node.success
    res <<- res
  }
  # get nodes across times
  nodes <- unique (unlist (llply (.data = clade.performance,
                                  .fun = function (x) as.vector(x$node))))
  # build res dataframe by adding first results
  res <- data.frame (.getNode (nodes[1]))
  colnames (res) <- nodes[1]
  nodes <- data.frame (node = nodes[-1])
  # add to res
  m_ply (.data = nodes, .fun = .addNode)
  rownames (res) <- time.intervals
  res
}

#' @name calcCladeStats
#' @title Calculate clade stats
#' @description Return a dataframe of stats for dataframe produced from getCladeSuccess().
#' @details The function calculates the following statistics for each clade:
#' 
#' * tot.size: the total number of tips in the clade through time
#' * max.size: the maximum number of tips in the clade in one time step
#' * time.span: the number of time steps the clade existed
#' * start: the starting time step
#' * end: the ending time step
#' * cm: the centre of mass of the clade
#' * cg: the centre of gravity of the clade
#' @template base_template
#' @param clades
#' @export
#' @examples
#' # example.var <- exampleFun (test.data)
#' @references
#' Gould, S. J., Gilinsky, N. L., & German, R. Z. (1987). Asymmetry of
#' lineages and the direction of evolutionary time. Science (New York, N.Y.),
#' 236(4807), 1437–1441.

calcCladeStats <- function (clades) {
  # Take dataframe of clade sizes through time, return stats
  #  for each clade
  .calc <- function (name) {
    clade <- clades[ ,as.character (name)]
    # basic stats
    tot.size <- sum (clade)
    max.size <- max (clade)
    if (tot.size == 0) {
      start <- end <- time.span <- NA
    } else {
      # get start and end
      start <- which (clade > 0)[1]
      end <- which (clade > 0)[sum (clade > 0)]
      # only count non-zero times
      clade <- clade[clade != 0]
      time.span <- length (clade)
    }
    if (time.span < 2 || is.na (time.span)) {
      cm <- cg <- NA
    } else {
      # calculate nominal times
      times <- seq (0, 1, 1/(length (clade) - 1))
      # calculate CM
      cm <- sum (clade * times)/sum (clade)
      # calculate CG
      cg <- sum (clade * (cm - times)^2)/sum (clade)
    }
    data.frame (tot.size, max.size, time.span, start, end, cm, cg)
  }
  mdply (.data = data.frame (name = colnames (clades)),
         .fun = .calc)
}

#' @name plotClades
#' @title Plot clade success through time
#' @description Return a geom_object of clade success through time.
#' @details This function will produce a geom_object() of the top k
#' clades by size through time. It takes a dataframe of clade successes
#' through time produced with getCladeSuccess().
#' 
#' By default the function will select only the biggest k clades by total
#' size. Alternatively, the user can select the clades of interest by
#' specifying their column indexes in the clades dataframe.
#' 
#' The user can calculate clade stats inside or outside the function.
#' @template base_template
#' @param clades
#' @param k number of clades to plot
#' @param clade.stats stats calculated using calcCladeStats (optional)
#' @param i indexes of user selected of clades to plot (optional)
#' @export
#' @examples
#' # example.var <- exampleFun (test.data)

plotClades <- function (clades, k=3, clade.stats=NULL, i=NULL) {
  # Return geom_object of k clade success through time
  .combine <- function (i) {
    # combine into single dataframe for ggplot
    n <- as.vector (clades[ ,i])
    t <- 1:length (n)
    c <- colnames (clades)[i]
    data.frame (n, t, c)
  }
  if (is.null (i)) {
    # If not i given, use k
    if (is.null (clade.stats)) {
      clade.stats <- calcCladeStats (clades)
    }
    # get top k clades by total size
    i <- order (clade.stats$tot.size, decreasing=TRUE)[1:k]
    # get their position in clades
    i <- which (clade.stats$name[i] %in% colnames (clades))
  }
  # combine clades into ggplot plot dataframe
  p.data <- mdply (.data=data.frame(i=i), .fun=.combine)[ ,-1]
  p <- ggplot (p.data, aes (x=t, y=n, colour=c)) +
    geom_line() + theme_bw() + theme (legend.position="none") +
    xlab ('Time') + ylab ('N') + ggtitle (paste0 ('K = ', k))
  return (p)
} 