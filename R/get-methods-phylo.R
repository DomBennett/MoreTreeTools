# All hidden functions for get-methods-visible.R using phylo class
.getTreeStats__phylo <- function (tree) {
  .get <- function (node, res) {
    # assign -- present
    next.edges <- which (tree$edge[ ,1] == node)
    next.nodes <- tree$edge[tree$edge[ , 1] == node, 2]
    if (length (next.nodes) == 0) {
      next.nodes <- next.edges <- NULL
    }
    res[[node]][['next.edges']] <- next.edges
    res[[node]][['next.nodes']] <- next.nodes
    res[[node]][['prev.edges']] <- which (tree$edge[ ,2] == node)
    res[[node]][['prev.nodes']] <- tree$edge[tree$edge[ , 2] == node, 1]
    # assign -- past
    res[[node]][['ascend.nodes']] <- res[[node]][['ascend.nodes']]
    res[[node]][['ascend.edges']] <- c (res[[node]][['ascend.edges']],
                                        res[[node]][['prev.edges']])
    # assign -- future
    if (is.null (next.nodes)) {
      # set for a tip
      descend.edges <- descend.nodes <- NULL
      children <- tree$tip.label[node]
      age <- 0
      pd <- 0
    } else {
      descend.edges <- next.edges
      descend.nodes <- next.nodes
      children <- NULL
      pd <- sum (tree$edge.length[next.edges])
      # setup and run for the next nodes
      for (next.node in next.nodes) {
        res[[next.node]] <- list ()
        res[[next.node]][['ascend.nodes']] <- c (res[[node]][['ascend.nodes']], node)
        res[[next.node]][['ascend.edges']] <- res[[node]][['ascend.edges']]
        res <- .get (next.node, res)
        descend.edges <- c (descend.edges,
                            res[[next.node]][['descend.edges']])
        descend.nodes <- c (descend.nodes,
                            res[[next.node]][['descend.nodes']])
        children <- c (children, res[[next.node]][['children']])
        pd <- pd + res[[next.node]][['pd']]
      }
      # only use one of the next nodes for age
      age <- root.age - sum (tree$edge.length [res[[node]][['ascend.edges']]])
    }
    res[[node]][['descend.edges']] <- descend.edges
    res[[node]][['descend.nodes']] <- descend.nodes
    res[[node]][['children']] <- children
    res[[node]][['n.children']] <- length (children)
    res[[node]][['age']] <- age
    res[[node]][['pd']] <- pd
    # return
    res
  }
  if (is.null (tree$edge.length)) {
    tree$edge.length <- rep (1, nrow (tree$edge))
  }
  # get age
  root.age <- .getSize__phylo (tree, 'rtt')
  # start with root node
  root.node <- length (tree$tip.label) + 1
  # init res list
  res <- list ()
  res[[root.node]] <- list ()
  res[[root.node]][['ascend.nodes']] <-
    res[[root.node]][['ascend.edges']] <- NULL
  .get (root.node, res)
}
.getNodeStats__phylo <- function (tree, nodes='all', ignore.tips=NULL) {
  .calc <- function (node) {
    children <- .getChildren__phylo (tree, node)
    children <- children[!children %in% ignore.tips]
    n.children <- length (children)
    age <- .getAge__phylo (tree, node=node)[ ,2]
    previous.node <- tree$edge[node == tree$edge[ ,2], 1]
    if (length (previous.node) == 0) {
      previous.node <- NA
    }
    data.frame (n.children, age, previous.node)
  }
  if (nodes[1] == 'all') {
    nodes <- 1:tree$Nnode + length (tree$tip.label)
  }
  loop.data <- data.frame (node=nodes, stringsAsFactors=FALSE)
  res <- mdply (.data=loop.data, .fun=.calc)
  if (!is.null (tree$node.label)) {
    res$node.label <- tree$node.label[nodes]
  }
  return (res)
}
.getEdgeStats__phylo <- function (tree, edges='all', ignore.tips=NULL) {
  .calc <- function (edge) {
    node.1 <- tree$edge[edge, 1]
    node.2 <- tree$edge[edge, 2]
    age.1 <- .getAge__phylo (tree, node=node.1)[ ,2]
    age.2 <- .getAge__phylo (tree, node=node.2)[ ,2]
    children <- .getChildren__phylo (tree, node.2)
    children <- children[!children %in% ignore.tips]
    n.children <- length (children)
    res <- data.frame (node.1, node.2, age.1, age.2,
                       n.children)
    return (res)
  }
  if (edges[1] == 'all') {
    edges <- 1:nrow (tree$edge)
  }
  loop.data <- data.frame (edge=edges, stringsAsFactors=FALSE)
  res <- mdply (.data=loop.data, .fun=.calc)
  if (!is.null (tree$edge.label)) {
    res$edge.label <- tree$edge.label[edges]
  }
  return (res)
}
.getOutgroup__phylo <- function (tree, tips) {
  # find the outgroup of a set of tips
  shrunk <- drop.tip (tree,
                      tip=tree$tip.label[!tree$tip.label %in% tips])
  shrunk$edge.length <- rep (1, nrow (shrunk$edge))
  dists <- colSums (cophenetic.phylo (shrunk))
  names (dists)[dists == max (dists)]
}
.getExtant__phylo <- function (tree, tol) {
  tip.ages <- .getAge__phylo (tree, 1:length (tree$tip.label))
  pull <- tip.ages[, 2] < tol
  tree$tip.label[pull]
}
.getChildren__phylo <- function(tree, node, display=FALSE) {
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
.getSize__phylo <- function (tree, type = c ('ntips', 'pd', 'rtt')) {
  type = match.arg (type)
  if (type == 'ntips') {
    return (length (tree$tip.label))
  } else if (type == 'pd') {
    return (sum (tree$edge.length))
  } else {
    age <- max (diag (vcv.phylo (tree)))
    if (!is.null (tree$root.edge)) {
      age <- age + tree$root.edge
    }
    return (age)
  }
}
.getAge__phylo <- function (tree, node='all', edge=NULL) {
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
    if (!is.numeric (node) || node > length (tree$tip.label) + tree$Nnode) {
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
.getSister__phylo <- function (tree, node = 'all') {
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
    nodes <- (length (tree$tip.label) + 2):(length (tree$tip.label) + tree$Nnode)
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
.getParent__phylo <- function (tree, node=NULL, tips=NULL, edges=NULL) {
  if (!is.null (node) & length (node) == 1) {
    if (!is.numeric (node)) {
      stop ('Node must be numeric')
    }
    if (node > length (tree$tip.label) + tree$Nnode) {
      stop ('Node not in tree')
    }
    if ((node == length (tree$tip.label) + 1) & is.rooted (tree)) {
      # if node is root, return it
      return (node)
    }
    return (tree$edge[tree$edge[ ,2] == node, 1])
  } else if (!is.null (tips)) {
    if (is.character (tips)) {
      # if tips are labels
      edges <- match (match (tips, tree$tip.label), tree$edge[,2])
    } else {
      # ... else they're numbers
      edges <- match (tips, tree$edge[,2])
    }
  } else if (!is.null (node)) {
    edges <- which (tree$edge[ ,2] %in% node)
  } else if (!is.null (edges)) {
    if (is.character (edges) & !is.null (tree$edge.label)) {
      # assume they are labels
      edges <- match (edges, tree$edge.label)
    }
  } else {
    stop ('Must provide either edges, tips or nodes argument')
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
.getEdges__phylo <- function (tree, node = NULL, tips = NULL, type = 1) {
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
.getNodes__phylo <- function (tree, node) {
  ## Find all nodes from given node to root
  base.node <- length (tree$tip.label) + 1
  nodes <- c ()
  while (node != base.node) {
    node <- tree$edge[tree$edge[ ,2] == node,1]
    nodes <- c (nodes, node)
  }
  nodes
}
.getClades__phylo <- function (tree) {
  if (is.null (tree$all.node)) {
    nodes <- 1:(tree$Nnode + length (tree$tip.label))
  } else {
    nodes <- tree$all.node
  }
  clades <- mlply (.data = data.frame (node = nodes),
                   .fun = .getChildren__phylo, tree)
  if (!is.null (tree$all.node.label)) {
    names (clades) <- tree$all.node.label
  }
  sizes <- ldply (.data = clades, .fun = length)[ ,2]
  clades <- clades[order (sizes, decreasing = TRUE)]
  nodes <- nodes[order (sizes, decreasing = TRUE)]
  list (clade.children = clades, clade.node = nodes)
}
.getNodeLabels__phylo <- function (tree, all = FALSE, datasource = 4) {
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
    children <- .getChildren__phylo (tree, node = i)
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
.getSubtrees__phylo <- function (tree, min.n, max.n, verbose = FALSE) {
  countChildren <- function (node) {
    # add node information to lists children and n
    these.children <- .getChildren__phylo (tree, node)
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
  ntips <- length (tree$tip.label)
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