#' @name getOutgroup
#' @title Find outgroup in a tree from tips
#' @description Return the outgroup tip of a tree from a vector of tips given.
#' @details Note, if polytomies occur in the tree, outgroup maybe a vector.
#' @template base_template
#' @param tips vector of tip labels
#' @export
#' @examples
#' # example.var <- exampleFun (test.data)

getOutgroup <- function (tree, tips) {
  # find the outgroup of a set of tips
  shrunk <- drop.tip (tree,
                      tip=tree$tip.label[!tree$tip.label %in% tips])
  shrunk$edge.length <- rep (1, nrow (shrunk$edge))
  dists <- colSums (cophenetic.phylo (shrunk))
  names (dists)[dists == max (dists)]
}

#' @name getExtant
#' @title Return extant tips
#' @description Return all tip labels whose tip age is 0.
#' @details No details
#' @template base_template
#' @param tol tolerance, below this is considered 0
#' @export
#' @examples
#' # example.var <- exampleFun (test.data)

getExtant <- function (tree, tol=1e-08) {
  # Return all extant tip labels
  tip.ages <- getAge (tree, 1:getSize (tree))
  pull <- tip.ages[, 2] < tol
  tree$tip.label[pull]
}

#' @name getChildren
#' @title Return descendant species from a node
#' @description Return all tip labels that descend from a specifed node.
#' @details No details
#' @template base_template
#' @param node number indicating node
#' @param display if true, a tree will be plotted with specifed node highlighted
#' @export
#' @examples
#' # example.var <- exampleFun (test.data)

getChildren <- function(tree, node, display = FALSE) {
  if (!is.numeric (node)) {
    stop("Node is not numeric!")
  }
  if (node > tree$Nnode + length (tree$tip.label)) {
    stop("Node is greater than the number of nodes in tree!")
  }
  if (node <= length (tree$tip.label)) {
    term.nodes <- node
  } else {
    term.nodes <- vector ()
    temp.nodes <- node
    while (length (temp.nodes) > 0) {
      connecting.nodes <- tree$edge[tree$edge[,1] %in% temp.nodes, 2]
      term.nodes <- c(term.nodes, connecting.nodes[connecting.nodes <=
                                                     length (tree$tip.label)])
      temp.nodes <- connecting.nodes[connecting.nodes > length(tree$tip.label)]
    }
  }
  children <- tree$tip.label[term.nodes]
  if (display) {
    tip.cols <- ifelse(tree$tip.label %in% children, "black", "grey")
    plot.phylo(tree, tip.color = tip.cols, show.tip.label = TRUE)
    nodelabels("node", node)
  }
  return (children)  
}

#' @name getSize
#' @title Return the size of a phylogenetic tree
#' @description Return either the number of tips, the total branch length
#' or the root to tip distance of a phylogenetic tree
#' @details I often need to quickly work out the size of a tree. This often means
#' writing \code{length (tree$tip.label)} over and over... but with this function
#' you can simply write \code{getSize(tree)}, that's a saving of 10 characters
#' (no sniffing matter for the hard core programmer) plus it looks neater. Not
#' to mention saving your brain time too.
#' @template base_template
#' @param type type of size to be calcualted: \code{ntips} (default), total branch length
#' \code{pd} or root to tip distance \code{rtt}
#' @export

getSize <- function (tree, type = c ('ntips', 'pd', 'rtt')) {
  type = match.arg (type)
  if (type == 'ntips') {
    return (length (tree$tip.label))
  } else if (type == 'pd') {
    return (sum (tree$edge.length))
  } else {
    return (max (diag (vcv.phylo (tree))))
  }
}

#' @name getAge
#' @title Return age of node or edge
#' @description Return data.frame of age of tip, internal node or age range of
#' an edge (so long as the tree provided is time calibrated with branch lengths).
#' @details First calculates the root to tip distance, then calculates
#' node age by subtracting this distance from the root to node distance.
#' If tree provided is not ultrametric, root to tip distance is the distance
#' from the root to the most distal tip.
#' If edge provided returns the max and min age of each node of edge.
#' If \code{node} equals 'all', will return a dataframe for all nodes.
#' @template base_template
#' @param node number(s) indicating node(s) (default 'all')
#' @param edge number(s) indicating edge(s) (default NULL)
#' @export
#' @examples
#' data ('catarrhines')
#' # when did humans and chimps split?
#' node <- getParent (catarrhines, tips=c ('Homo sapiens', 'Pan troglodytes'))
#' (getAge (catarrhines, node=node))  # 9.7MYA according to this tree

getAge <- function (tree, node='all', edge=NULL) {
  run <- function (node) {
    # function for calculating age of node
    term.node <- length (tree$tip.label) + 1
    # if it's the root node, return tree.age
    if (term.node == node) {
      return (tree.age)
    }
    # if it's a tip return 0
    if (node < length (tree$tip.label) & is.ultrametric (tree)) {
      return (0)
    }
    # else find all its edge.lengths and subtract from tree.age
    edges <- c ()
    while (node != term.node) {
      edges <- c (edges, which (tree$edge[ ,2] == node))
      node <- tree$edge[tree$edge[ ,2] == node, 1]
    }
    return (tree.age - sum (tree$edge.length[edges]))
  }
  runEdge <- function (edge) {
    max.age <- run (node=tree$edge[edge, 1])
    min.age <- run (node=tree$edge[edge, 2])
    data.frame (max.age, min.age)
  }
  # SAFETY CHECK
  if (node[1] != 'all') {
    if (!is.numeric (node) || node > getSize (tree) + tree$Nnode) {
      stop ('Node must either be `all` or a number specifying a node in tree.')
    }
  }
  if (!is.null (edge)) {
    if (!is.numeric (edge) || edge > nrow (tree$edge)) {
      stop ('Edge must be a numeric specifying an edge row in tree$edge')
    }
  }
  tree.age <- max (diag (vcv.phylo (tree)))
  if (!is.null (edge)) {
    res <- mdply (.data = data.frame (edge=edge), .fun = runEdge)
  } else if (node[1] != 'all') {
    res <- mdply (.data = data.frame (node=node), .fun = run)
    colnames (res)[2] <- 'age'
  } else {
    # else node == all, run on all nodes
    nodes <- 1:(length (tree$tip.label) + tree$Nnode)
    res <- mdply (.data = data.frame (node = nodes), .fun = run)
    colnames (res)[2] <- 'age'
  }
  return (res)
}

#' @name getSister
#' @title Return sister node
#' @description Return sister node of node(s) given, default will
#' return sister nodes of all internal nodes.
#' @details A sister node is defined as the other node descending
#' from the parent node. This functions finds all sister nodes.
#' It does not handle polytomies.
#' @template base_template
#' @param node number indicating node
#' @export
#' @examples
#' # example.var <- exampleFun (test.data)

getSister <- function (tree, node = 'all') {
  .getSister <- function (node) {
    inner.node <- tree$edge[match(node, tree$edge[,2]),1]
    # find all nodes and edges descending from inner node
    d.edges <- which(tree$edge[,1] %in% inner.node)
    d.nodes <- tree$edge[d.edges, 2]
    # remove starting node
    nearest.node <- d.nodes[!d.nodes %in% node][1]
    nearest.node
  }
  if (node[1] == 'all') {
    nodes <- (getSize (tree) + 2):(getSize (tree) + tree$Nnode)
    res <- mdply (.data = data.frame (node = nodes),
                  .fun = .getSister)
    colnames (res)[2] <- 'sister'
  } else if (length (node) > 1) {
    res <- mdply (.data = data.frame (node = node),
                  .fun = .getSister)
    colnames (res)[2] <- 'sister'
  } else {
    res <- .getSister (node)
  }
  return (res)
}

#' @name getParent
#' @title Return parent node
#' @description Return parent node of a node or a vector of tips
#' @template base_template
#' @param node number indicating node
#' @param tips tip labels or numbers

getParent <- function (tree, node=NULL, tips=NULL) {
  if (is.null (tips) & is.null (node)) {
    stop ('Must provide either tips or nodes argument')
  }
  if (!is.null (node)) {
    if (!is.numeric (node)) {
      stop ('Node must be numeric')
    }
    if (node > getSize (tree) + tree$Nnode) {
      stop ('Node not in tree')
    }
    if ((node == getSize (tree) + 1) & is.rooted (tree)) {
      # if node is root, return it
      return (node)
    }
    return (tree$edge[tree$edge[ ,2] == node, 1])
  }
  if (is.character (tips)) {
    # if tips are labels
    edges <- match (match (tips, tree$tip.label), tree$edge[,2])
  } else {
    # ... else they're numbers
    edges <- match (tips, tree$edge[,2])
  }
  end.nodes <- tree$edge[edges, 1]
  term.node <- length (tree$tip.label) + 1
  while (TRUE){
    if (sum (end.nodes[1] == end.nodes) == length (end.nodes)){
      break
    }
    end.nodes <- sort (end.nodes, TRUE)
    start.node <- end.nodes[1]
    edge <- match (start.node, tree$edge[,2])
    end.node <- tree$edge[edge,1]
    edges <- c(edges, edge)
    end.nodes <- c(end.nodes[!end.nodes %in% start.node], end.node)
  }
  return (end.nodes[1])
}

#' @name getEdges
#' @title Return edges
#' @description Return all children edges from a node,
#'  or all edges connected to tips based on three different methods
#' @details If node specified, all edges that descend from specified node are returned
#' 
#' Otherwise, if tips are specified edges are returned based on one of three methods:
#' Type:
#' 
#' 1 -- phylogeny consisting solely of the tips, default
#' 
#' 2 -- edges from taxon tips to terminal node
#' 
#' 3 -- edges unique to tips
#' @template base_template
#' @param node number indicating node
#' @param tips tip labels or numbers
#' @param type 1, 2 or 3 indicating which method to use
#' @export
#' @examples
#' # example.var <- exampleFun (test.data)

getEdges <- function (tree, node = NULL, tips = NULL, type = 1) {
  if (!is.null (node)) {
    ## Find all edges from given node to tips
    edges <- c ()
    while (TRUE) {
      bool <- tree$edge[ ,1] %in% node
      if (sum (bool) > 1) {
        node <- tree$edge[bool, 2]
        edges <- c (edges, which (bool))
      } else {
        break
      }
    }
    return (edges)
  }
  # otherwise, use tips...
  if (!type %in% c (1,2,3)) {
    stop ("Type must be an integer: 1, 2 or 3.")
  }
  if (length (tips) == length (tree$tip.label)){
    return (tree$edge)
  }
  if (type == 1 & length (tips) == 1){
    stop ("length (tips) == 1 :
         Cannot return a single edge for type 1.")
  }
  # start at the tips and step back into the phylogeny ...
  # ...add all connecting edges to a vector...
  # stop when all paths have met at the same node (type = 1)
  # or when all paths have reached the root node (type = 2)
  # or when all the nodes are unique (type = 3)
  if (is.character (tips)) {
    # if tips are labels
    edges <- match (match (tips, tree$tip.label), tree$edge[,2])
  } else {
    # ... else they're numbers
    edges <- match (tips, tree$edge[,2])
  }
  end.nodes <- tree$edge[edges, 1]
  term.node <- length (tree$tip.label) + 1
  if (all (end.nodes %in% term.node)) {
    return (edges)
  } else {
    if (type == 3){
      while (any (duplicated (end.nodes))){
        start.node <- end.nodes[duplicated(end.nodes)][1]
        if (sum (tree$edge[,1] %in% start.node) == sum (end.nodes %in% start.node)){
          edge <- match (start.node, tree$edge[,2])
          end.node <- tree$edge[edge,1]
          edges <- c(edges, edge)
          end.nodes <- c(end.nodes[!end.nodes %in% start.node], end.node)
        } else {
          end.nodes <- end.nodes[end.nodes != start.node]
        }
      }
    } else {
      while (TRUE){
        if (type == 2){
          if (sum (term.node == end.nodes) == length (end.nodes)){
            break
          }
        } else {
          if (sum (end.nodes[1] == end.nodes) == length (end.nodes)){
            break
          }
        }
        end.nodes <- sort (end.nodes, TRUE)
        start.node <- end.nodes[1]
        edge <- match (start.node, tree$edge[,2])
        end.node <- tree$edge[edge,1]
        edges <- c(edges, edge)
        end.nodes <- c(end.nodes[!end.nodes %in% start.node], end.node)
      }
    }
    return (edges)
  }
}

#' @name getNodes
#' @title Return parent nodes
#' @description Return all nodes that connect the specified node to the root
#' @details No details
#' @template base_template
#' @param node number indicating node
#' @export
#' @examples
#' # example.var <- exampleFun (test.data)

getNodes <- function (tree, node) {
  ## Find all nodes from given node to root
  base.node <- length (tree$tip.label) + 1
  nodes <- c ()
  while (node != base.node) {
    node <- tree$edge[tree$edge[ ,2] == node,1]
    nodes <- c (nodes, node)
  }
  nodes
}

#' @name getClades
#' @title Return all children of each node
#' @description All children for each node in the tree are returned.
#' @details Returns a list of lists: 'clade.children' contains all the children of each clade,
#'  'clade.node' contains the clade number of each clade in clade.children.
#' @template base_template
#' @export
#' @examples
#' # example.var <- exampleFun (test.data)

getClades <- function (tree) {
  if (is.null (tree$all.node)) {
    nodes <- 1:(tree$Nnode + length (tree$tip.label))
  } else {
    nodes <- tree$all.node
  }
  clades <- mlply (.data = data.frame (node = nodes),
                   .fun = getChildren, tree)
  if (!is.null (tree$all.node.label)) {
    names (clades) <- tree$all.node.label
  }
  sizes <- ldply (.data = clades, .fun = length)[ ,2]
  clades <- clades[order (sizes, decreasing = TRUE)]
  nodes <- nodes[order (sizes, decreasing = TRUE)]
  list (clade.children = clades, clade.node = nodes)
}

#' @name getNodeLabels
#' @title Return node labels based on online taxonomic database
#' @description Return names of each node in tree based on searching tip labels
#'  in the Global Names Resolver \url{http://resolver.globalnames.org/}
#' @details For each node, all the children are searched, the taxonomic lineages returned and
#' then searched to find the lowest shared name.
#' All the tip labels are searched against a specified taxonomic database through the GNR.
#'  The default is NCBI (4). For a list of other possible datasources see:
#'   \url{http://resolver.globalnames.biodinfo.org/data_sources}
#' @template base_template
#' @param all count tips as nodes, default False
#' @param datasource a number indicating the GNR datasource to search against
#' @export
#' @examples
#' # example.var <- exampleFun (test.data)

getNodeLabels <- function (tree, all = FALSE, datasource = 4) {
  # Use GNR to label all nodes in a phylogeny
  # first replace all _ with spaces
  tree$tip.label <- gsub ('_', ' ', tree$tip.label)
  taxa.res <- taxaResolve (tree$tip.label, datasource = datasource)
  nodes <- 1:(length (tree$tip.label) + tree$Nnode)
  node.label <- rep (NA, length (nodes))
  # for tips use the first word of the name
  node.label[1:length (tree$tip.label)] <-
    unlist (lapply (strsplit (tree$tip.label, "\\s+"), function (x) x[1]))
  for (i in (length (tree$tip.label) + 1):length (nodes)) {
    children <- getChildren (tree, node = i)
    genus.names <- unlist (lapply (strsplit (children, " "), function (x) x[1]))
    if (all (genus.names == genus.names[1])) {
      node.label[i] <- genus.names[1]
    } else {
      lineages <- as.character (taxa.res[taxa.res$search.name %in% children,
                                         "lineage"])
      lineages <- strsplit (lineages, "\\|")
      lineages <- lineages[!is.na (lineages)]
      if (length (lineages) > 0) {
        node.label[i] <- .findClade (lineages)
      }
    }
  }
  if (all) {
    return (node.label)
  } else {
    return (node.label[-(1:length(tree$tip.label))])
  }
}

#' @name getSubtrees
#' @title Return non-redundant trees of clades from a phylogenetic tree
#' @description Extract subtrees from a phylogenetic tree within a specified range of
#' number of tips, where each subtree returned has tips unique to it.
#' @details This functions searches through all the nodes in a tree, it
#' identifies all nodes that have numbers of descendents in the range
#' specified by the user, it then returns the largest nodes for which all tips
#' are unique to each subtree.
#' @template base_template
#' @param min.n the minimum number of tips to return for each clade
#' @param max.n the maximum number of tips to return for each clade
#' @param verbose true or false
#' @export
#' @examples
#' # example.var <- exampleFun (test.data)

getSubtrees <- function (tree, min.n, max.n, verbose = FALSE) {
  countChildren <- function (node) {
    # add node information to lists children and n
    these.children <- getChildren (tree, node)
    this.n <- length (these.children)
    if (this.n <= max.n & this.n >= min.n) {
      children <<- c (children, list (these.children))
      n <<- c (n, this.n)
      node.number <<- c (node.number, node)
    }
  }
  checkNode <- function (i, these.names) {
    # return True if names do not overlap
    those.names <- children[[i]]
    !any (these.names %in% those.names)
  }
  # how many nodes to loop through
  ntips <- getSize (tree)
  nodes <- (ntips + 1):(ntips + tree$Nnode)
  # quick return if it's already right size
  if (ntips <= max.n && ntips >= min.n) {
    if (verbose) {
      cat (paste0 ('\nTree already within min and max n'))
    }
    tree <- list (tree)
    class (tree) <- 'multiPhylo'
    return (tree)
  }
  # create a list for tip names for each clade
  children <- list ()
  # create a vector of n tips and node number for each clade
  n <- node.number <- NULL
  # loop through nodes writing info to children and n
  m_ply (.data = data.frame (node = nodes), .fun = countChildren)
  if (is.null (n)) {
    if (verbose) {
      cat (paste0 ('\nNo subtreess found between [', min.n,
                   '] and [', max.n,']'))
    }
    return (NULL)
  }
  if (length (children) == 1) {
    # if length is 1, then only one clade matched
    tree <- extract.clade (tree, node = node.number)
    tree <- list (tree)
    class (tree) <- 'multiPhylo'
    return (tree)
  }
  # out of those clades, find a non-redundant combination
  this <- 1
  while (this <= length (n)) {
    # work out overlap between this node and other nodes
    those <- (1:length (n))[-this]
    these.names <- children[[this]]
    bool <- mdply (.data = data.frame (i = those), .fun = checkNode,
                   these.names)[ ,2]
    # if this node's n is greater than its overlapping 'those nodes' keep,
    #  else drop
    this.n <- n[this]
    these.n <- n[those[!bool]]
    if (any (this.n <= these.n)) {
      n <- n[-this]
      children <- children[-this]
      node.number <- node.number[-this]
    } else {
      this <- this + 1
    }
  }
  # extract clades and return as trees
  trees <- list ()
  for (each in node.number) {
    clade.tree <- extract.clade (tree, node = each)
    trees <- c (trees, list (clade.tree))
  }
  # return as multiPhylo
  class (trees) <- 'multiPhylo'
  return (trees)
}