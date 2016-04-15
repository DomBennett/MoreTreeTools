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
#' Setting 'ind' (or independent) as TRUE will implement a counting
#' method that attempts to count descendents unique to the clade by
#' preventing double counting in decsendent clades. Independent
#' counting works by limiting the maximum clade count to the minimum
#' of the descendant counts of the two daughter clades. If one of the
#' daughter clades continues to speciate while the other goes extinct
#' then the parent clade of the two will be limited to 0 rather than
#' the size of a single daughter clade.
#' 
#' 'time.intervals' is a vector of time steps represented by each tree.
#' Default is 1 to length of number of trees.
#' 
#' Setting 'count.extincts' to TRUE will count both extinct and extant
#' tips through time.
#' 
#' @template base_template
#' @param trees
#' @param ind boolean, use independent counting? Deafult FALSE
#' @param count.extincts boolean, count extinct and extant tips? Default FALSE
#' @param time.intervals vector of time represented by each tree in trees
#' @export
#' @examples
#' library (MoreTreeTools)
#' # simulate an EDBMM tree (this will take time!)
#' trees <- runEDBMM (birth=2, death=1, stop.at=200, record=TRUE,
#'                    progress.bar='time')
#' # get raw number of descendents by node for trees
#' clades <- getCladeSuccess (trees)
#' # get independent number of descendents by node
#' clades.ind <- getCladeSuccess (trees, ind=TRUE)
#' # calculate stats by clade
#' clade.stats <- calcCladeStats (clades)
#' clade.ind.stats <- calcCladeStats (clades.ind)
#' # plot the trajectories of the 5 biggest clades
#' plotClades (clades=clades, clade.stats=clade.stats, k=5)
#' # plot clades through time for all clades that went extinct
#' # by the end of the simulation (there are 11 time steps in the simulation)
#' selection <- which (clade.stats$end < 11 &
#'                       clade.stats$max.size > 2)
#' plotClades (clades=clades, i=selection)
#' plotClades (clades=clades, i=selection, merge=TRUE)
#' # compare with ind clades
#' selection <- which (clade.ind.stats$end < 11 &
#'                       clade.ind.stats$max.size > 2)
#' plotClades (clades=clades.ind, i=selection)
#' plotClades (clades=clades.ind, i=selection, merge=TRUE)
#' # it is difficult to detect a difference in clade shape
#' # with so few time steps. Instead let's load some generated
#' # test clades
#' data ('testclades')
#' # this loads 'unif.clades' and 'arch.clades'
#' plotClades (clades=unif.clades, i=1:ncol (unif.clades), merge=TRUE)
#' plotClades (clades=arch.clades, i=1:ncol (arch.clades), merge=TRUE)
#' # now we can see a difference

getCladeSuccess <- function (trees, ind=FALSE, count.extincts=FALSE,
                             time.intervals=1:length(trees)) {
  # Return dataframe of clades through time from list of trees
  # with named nodes
  .get <- function (i) {
    .countChildren (trees[[i]], ind, count.extincts)
  }
  clade.successes <- mlply (.data = data.frame (i = 1:length (trees)),
                            .fun = .get)
  .reformat (clade.successes, time.intervals)
}
.countChildren <- function (tree, ind, count.extincts) {
  # Count the number of extant children for every node
  .count <- function (node.label=NULL, node=NULL) {
    if (!is.null (node.label)) {
      node <- which (tree$node.label == node.label)
      node <- length (tree$tip.label) + node
    }
    node.children <- getChildren (tree, node)
    extant <- node.children[node.children %in% extants]
    return (length (extant))
  }
  .countind <- function (node.label) {
    node <- which (tree$node.label == node.label)
    node <- length (tree$tip.label) + node
    d.nodes <- tree$edge[which (tree$edge[ ,1] == node),2]
    d1 <- .count (node=d.nodes[1])
    d2 <- .count (node=d.nodes[2])
    # return weighted clade number
    d1 + d2 - abs (d1 - d2)
  }
  # all internal nodes
  node.labels <- tree$node.label
  if (count.extincts) {
    extants <- tree$tip.label
  } else {
    # make sure only extant clades are counted
    extants <- getExtant (tree)
  }
  if (ind) {
    res <- mdply (.data = data.frame (node.label = node.labels),
                  .fun = .countind)
  } else {
    res <- mdply (.data = data.frame (node.label = node.labels),
                  .fun = .count)
  }
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
#' tot.size: the total number of tips in the clade through time
#' 
#' max.size: the maximum number of tips in the clade in one time step
#' 
#' time.span: the number of time steps the clade existed
#' 
#' start: the starting time step
#' 
#' end: the ending time step
#' 
#' cm: the centre of mass of the clade
#' 
#' cg: the centre of gyration of the clade
#' 
#' @param clades dataframe where columns are clades and rows are numbers
#'  of descendents by time step as returned by getCladeSuccess()
#' @export
#' @examples
#' library (MoreTreeTools)
#' # simulate an EDBMM tree (this will take time!)
#' trees <- runEDBMM (birth=2, death=1, stop.at=200, record=TRUE,
#'                    progress.bar='time')
#' # get raw number of descendents by node for trees
#' clades <- getCladeSuccess (trees)
#' # get independent number of descendents by node
#' clades.ind <- getCladeSuccess (trees, ind=TRUE)
#' # calculate stats by clade
#' clade.stats <- calcCladeStats (clades)
#' clade.ind.stats <- calcCladeStats (clades.ind)
#' # plot the trajectories of the 5 biggest clades
#' plotClades (clades=clades, clade.stats=clade.stats, k=5)
#' # plot clades through time for all clades that went extinct
#' # by the end of the simulation (there are 11 time steps in the simulation)
#' selection <- which (clade.stats$end < 11 &
#'                       clade.stats$max.size > 2)
#' plotClades (clades=clades, i=selection)
#' plotClades (clades=clades, i=selection, merge=TRUE)
#' # compare with ind clades
#' selection <- which (clade.ind.stats$end < 11 &
#'                       clade.ind.stats$max.size > 2)
#' plotClades (clades=clades.ind, i=selection)
#' plotClades (clades=clades.ind, i=selection, merge=TRUE)
#' # it is difficult to detect a difference in clade shape
#' # with so few time steps. Instead let's load some generated
#' # test clades
#' data ('testclades')
#' # this loads 'unif.clades' and 'arch.clades'
#' plotClades (clades=unif.clades, i=1:ncol (unif.clades), merge=TRUE)
#' plotClades (clades=arch.clades, i=1:ncol (arch.clades), merge=TRUE)
#' # now we can see a difference
#' @references
#' Gould, S. J., Gilinsky, N. L., & German, R. Z. (1987). Asymmetry of
#' lineages and the direction of evolutionary time. Science (New York, N.Y.),
#' 236(4807), 1437–1441.

calcCladeStats <- function (clades) {
  # Take dataframe of clade sizes through time, return stats
  #  for each clade
  .calc <- function (cid) {
    clade <- clades[ ,as.character (cid)]
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
  mdply (.data = data.frame (cid = colnames (clades),
                             stringsAsFactors=FALSE),
         .fun = .calc)
}

#' @name plotClades
#' @title Plot clade success through time
#' @description Return a geom_object of clade success through time.
#' @details This function will return a geom_object for plotting
#' the trajectories of clade success through time from a clades
#' dataframe generated by getCladeSuccess().
#' 
#' By default, the largest three by maximum size clades will be plotted.
#' This number can be changed with the N parameter. If using top N, for
#' dataframes of large sizes, it is computationally more efficient to
#' generate the clade.stats dataframe before calling this function with
#' calcCladeStats().
#' 
#' The user can also choose to plot the trajectories of specific clades
#' in the clades dataframe using 'cids' e.g. to show only clades that
#' have gone extinct or only clades of a certain size.
#' 
#' The user can also merge the clade trajectories by setting 'merge' as
#' TRUE. This will take the trajectory of each clade and normalise their
#' time and number of descendents by putting them onto the same scale.
#' Normalised time is calculated as the time when the clade appeared (0)
#' to when it went it extinct (1). Normalise number of descendents is
#' rescaled by setting the smallest size of the clade to 0 and the 
#' largest size to 1. Black line indicates the mean while red dashed
#' indicate the SD.
#' 
#' @param clades
#' @param N top number of clades to plot
#' @param clade.stats stats calculated using calcCladeStats (optional)
#' @param cids clade ids of user selected clades to plot (optional)
#' @param legend boolean, display legend? Default FALSE
#' @param merge boolean, plot normalised clade trajectory? Default FALSE
#' @export
#' @examples
#' library (MoreTreeTools)
#' # simulate an EDBMM tree (this will take time!)
#' trees <- runEDBMM (birth=2, death=1, stop.at=200, record=TRUE,
#'                    progress.bar='time')
#' # get raw number of descendents by node for trees
#' clades <- getCladeSuccess (trees)
#' # get independent number of descendents by node
#' clades.ind <- getCladeSuccess (trees, ind=TRUE)
#' # calculate stats by clade
#' clade.stats <- calcCladeStats (clades)
#' clade.ind.stats <- calcCladeStats (clades.ind)
#' # plot the trajectories of the 5 biggest clades
#' plotClades (clades=clades, clade.stats=clade.stats, k=5)
#' # plot clades through time for all clades that went extinct
#' # by the end of the simulation (there are 11 time steps in the simulation)
#' selection <- which (clade.stats$end < 11 &
#'                       clade.stats$max.size > 2)
#' plotClades (clades=clades, i=selection)
#' plotClades (clades=clades, i=selection, merge=TRUE)
#' # compare with ind clades
#' selection <- which (clade.ind.stats$end < 11 &
#'                       clade.ind.stats$max.size > 2)
#' plotClades (clades=clades.ind, i=selection)
#' plotClades (clades=clades.ind, i=selection, merge=TRUE)
#' # it is difficult to detect a difference in clade shape
#' # with so few time steps. Instead let's load some generated
#' # test clades
#' data ('testclades')
#' # this loads 'unif.clades' and 'arch.clades'
#' plotClades (clades=unif.clades, i=1:ncol (unif.clades), merge=TRUE)
#' plotClades (clades=arch.clades, i=1:ncol (arch.clades), merge=TRUE)
#' # now we can see a difference

plotClades <- function (clades, N=3, clade.stats=NULL, cids=NULL,
                        legend=FALSE, merge=FALSE) {
  # Internals
  .combine <- function (cid) {
    # combine into single dataframe for ggplot
    n <- as.vector (clades[ ,cid])
    t <- 1:length (n)
    data.frame (n, t, c=cid)
  }
  .normalise <- function (cid) {
    n <- as.vector (clades[ ,cid])
    n <- n[n != 0]
    n <- n / min (n)
    n <- n / max(n)
    t <- seq (0, 1, length.out=length (n))
    data.frame (n, t)
  }
  if (is.null (cids)) {
    if (merge) {
      # if no cids given and merge, use all clades
      cids <- colnames (clades)
    } else {
      # otherwise use top N clades
      if (is.null (clade.stats)) {
        clade.stats <- calcCladeStats (clades)
      }
      # get top N clades by total size
      cids <- order (clade.stats$tot.size, decreasing=TRUE)[1:N]
      # get their position in clades
      cids <- which (clade.stats$cid %in% colnames (clades))
    }
  }
  gt <- paste0 ('N = ', length (cids))
  p.data <- data.frame(cid=cids, stringsAsFactors=FALSE)  # avoid levels
  if (merge) {
    # if merged into single plot
    p.data <- mdply (.data=data.frame(p.data), .fun=.normalise)[ ,-1]
    dps <- round (log(x=nrow (clades), base=10))
    p.data$t <- signif (p.data$t, dps)
    p.data <- ddply (p.data, .(t), summarise,
                     mean.n=mean(n),
                     upper.sd=mean(n) + sd (n),
                     lower.sd=mean(n) - sd (n))
    p <- ggplot (p.data) +
      geom_line (aes (x=t, y=mean.n), colour='black') +
      geom_line (aes (x=t, y=upper.sd), colour='red', lty=3) +
      geom_line (aes (x=t, y=lower.sd), colour='red', lty=3) +
      theme_bw() + xlab ('Normalised time step') + ylab ('Normalised N')
  } else {
    # combine clades into ggplot plot dataframe
    p.data <- mdply (.data=p.data, .fun=.combine)[ ,-1]
    p <- ggplot (p.data, aes (x=t, y=n, colour=c)) +
      geom_line() + theme_bw() +
      xlab ('Time step') + ylab ('N')
  }
  if (!legend) {
    p <- p + theme (legend.position="none")
  }
  return (p + ggtitle(gt))
}