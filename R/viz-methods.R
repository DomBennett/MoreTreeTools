#' @name commplot
#' @title Plot a community or traits on a phylogenetic tree
#' @description Plot community or trait data on phylogeny to visualise
#' differences in structure between traits or sites. Use colours to distinguish
#' groups and alpha to distinguish abundances (works like rainbow())
#' @details No details
#' @template base_template
#' @param cmatrix community/trait data matrix (cols taxa, rows sites)
#' @param groups trait or site groups
#' @export
#' @examples
#' # Generate a balanced tree
#' tree <- compute.brlen (stree (16, 'balanced'))
#' # Generate random community data, 2 of 3 different communities
#' cmatrix <- matrix (round (runif (16*6, 0, 2)), nrow = 6, byrow = TRUE)
#' colnames (cmatrix) <- tree$tip.label
#' rownames (cmatrix) <- paste0 ('site_', 1:6)
#' # Plot community, specify group identitiy of each of the 6 sites
#' commplot (cmatrix, tree, groups = c (1,1,2,2,3,3))

commplot <- function(cmatrix, tree, groups = rep(1, nrow(cmatrix)),
                     no.margin = TRUE, ...){
  # ultrametricize tree FIRST
  tree <- compute.brlen(tree, method="Grafen")
  # plot phylogeny, allow space for points
  tree$edge.length <- tree$edge.length/getSize (tree, 'rtt')
  # make all trees the same length before plotting i.e. all branches from terminal
  # node to tip equal 1
  # for some weird reason the rules of plotting are dyanmic!
  if (nrow(cmatrix) < 20) {
    variable.max <- 1 + nrow(cmatrix)/20
    spacing.start <- 0.55
    spacing.i <- 0.05
  } else {
    variable.max <- nrow(cmatrix)/10
    spacing.i <- 0.1 - 1/nrow(cmatrix)
    spacing.start <- 0.5 + spacing.i
  }
  plot(tree, no.margin = no.margin, show.tip.label = FALSE,
       x.lim = c(0, variable.max), ...)
  # generate alphas based on abundances
  n <- length(unique(groups))
  hs <- seq.int(0, 1 + max(1, n - 1)/n, length.out = n)%%1
  alphas <- cmatrix/max(cmatrix)
  # loop init
  ntips <- length(tree$tip.label)
  spacing <- spacing.start
  group <- groups[1]
  j <- 1
  # loop through sites and plot points for present species
  for(i in 1:nrow(cmatrix)){
    j <- ifelse(group == groups[i], j, j + 1)
    pull <- as.logical(cmatrix[i,])
    taxa <- tree$tip.label[pull]
    abunds <- alphas[i, pull]
    tiplabels(tip = match(taxa, tree$tip.label),
              pch = 19, adj = spacing, col = hsv(rep(hs[j], ntips), 1, 1, abunds))
    spacing <- spacing + spacing.i
    group <- groups[i]
  }
}