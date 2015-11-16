#' @name getTreeStats
#' @title Return large list of trees stats by node
#' @description Returns a list object of nodes with information
#' on each node -- see details.
#' @details Many functions require the navigation of a tree e.g.
#' finding all descending nodes, counting number of children by node
#' etc. \code{MoreTreeTools} has many functions for performing these
#' tasks, but in cases where these function need to be run multiple
#' times it can be more efficient to generate this information beforehand.
#' 
#' A list object is returned, with stats for each node. Stats
#' generated are:
#' \itemize{
#'   \item All nodes from root to node, \code{ascend.nodes}
#'   \item All edges from root to edge, \code{ascend.edges}
#'   \item Direct previous node connecting to node, \code{prev.nodes}
#'   \item Direct previous edges connecting to node, \code{prev.edges}
#'   \item Tips descending from node, \code{children}
#'   \item Number of tips descending from node, \code{n.children}
#'   \item Age of node, \code{age}
#'   \item Phylogenetic diversity of node, \code{pd}
#' }
#' 
#' Note, this function assumes trees are rooted. If no branch lengths are
#' given, the function sets all branch lengths to 1. The function is
#' recursive, for very large trees (100,000 + tips) a stack overflow warning
#' may appear. Update \code{options(expressions=)} to change the stack limit.
#' @param tree
#' @export
#' @examples
#' # get hominoids tree
#' data ('hominoids')
#' # plot with labels
#' plot (hominoids);nodelabels();edgelabels();axisPhylo();
#' # generate tree.stats
#' tree.stats <- getTreeStats (hominoids)
#' # corroborate with plot
#' tree.stats[[19]][['prev.edges']]  # the previous edge for node 19 is 2
#' tree.stats[[20]][['children']]  # all the 'great' apes descend from node 20

getTreeStats <- function (tree) {
  .fun <- .phyloOrNodelist (tree, .getTreeStats__phylo,
                            .getTreeStats__NodeList)
  .fun (tree)
}

#' @name getNodeStats
#' @title Get node statistics
#' @description Return dataframe of number of children, age and
#' previous connecting node of every node in tree.
#' @details Assumes tree is rooted, bifurcating and time callibrated.
#' 
#' Use ignore.tips to specify tips of tree to ignore when counting
#' children.
#' @template base_template
#' @param nodes numeric vector of nodes, default 'all'
#' @param ignore.tips tips to ignore when counting, default NULL
#' @export
#' @examples
#' tree <- rtree (100)
#' node.stats <- getNodeStats (tree)

getNodeStats <- function (tree, nodes='all', ignore.tips=NULL) {
  .fun <- .phyloOrNodelist (tree, .getNodeStats__phylo,
                            .getNodeStats__NodeList)
  .fun (tree, nodes=nodes, ignore.tips=ignore.tips)
}

#' @name getEdgeStats
#' @title Get edge statistics
#' @description Return dataframe of first node of edge, second node of edge,
#' their ages and the number of children descending from the edge.
#' @details Assumes tree is rooted, bifurcating and time callibrated.
#' 
#' Use ignore.tips to specify tips of tree to ignore when counting
#' children.
#' @template base_template
#' @param edges numeric vector of edges, default 'all'
#' @param ignore.tips tips to ignore when counting, default NULL
#' @export
#' @examples
#' tree <- rtree (100)
#' edge.stats <- getEdgeStats (tree)

getEdgeStats <- function (tree, edges='all', ignore.tips=NULL) {
  .fun <- .phyloOrNodelist (tree, .getEdgeStats__phylo,
                            .getEdgeStats__NodeList)
  .fun (tree, edges=edges, ignore.tips=ignore.tips)
}

#' @name getOutgroup
#' @title Get outgroup in a tree from tips
#' @description Return the outgroup tip of a tree from a vector of tips given.
#' @details Note, if polytomies occur in the tree, outgroup maybe a vector.
#' @template base_template
#' @param tips vector of tip labels
#' @export
#' @examples
#' # example.var <- exampleFun (test.data)

getOutgroup <- function (tree, tips) {
  .fun <- .phyloOrNodelist (tree, .getOutgroup__phylo,
                            .getOutgroup__NodeList)
  .fun (tree, tips)
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
  .fun <- .phyloOrNodelist (tree, .getExtant__phylo,
                            .getExtant__NodeList)
  .fun (tree, tol)
}

#' @name getChildren
#' @title Get descendant species from a node
#' @description Return all tip labels that descend from a specifed node.
#' @details No details
#' @template base_template
#' @param node number indicating node
#' @param display if true, a tree will be plotted with specifed node highlighted (phylo only)
#' @export
#' @examples
#' # example.var <- exampleFun (test.data)

getChildren <- function (tree, node, display=FALSE) {
  # Return all child tips
  .fun <- .phyloOrNodelist (tree, .getChildren__phylo,
                            .getChildren__NodeList)
  .fun (tree, node, display=display)
}

#' @name getSize
#' @title Get the size of a phylogenetic tree
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
  .fun <- .phyloOrNodelist (tree, .getSize__phylo,
                            .getSize__NodeList)
  .fun (tree, type)
}

#' @name getAge
#' @title Get age of node or edge
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
  .fun <- .phyloOrNodelist (tree, .getAge__phylo,
                            .getAge__NodeList)
  .fun (tree, node=node, edge=edge)
}

#' @name getSister
#' @title Get sister node
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
  .fun <- .phyloOrNodelist (tree, .getSister__phylo,
                            .getSister__NodeList)
  .fun (tree, node)
}

#' @name getParent
#' @title Get parent node
#' @description Return parent node of nodes, edges or a vector of tips
#' @template base_template
#' @param node vector of number indicating node
#' @param tips vector of tip labels or numbers
#' @param edges vector of edge labels or edge numbers

getParent <- function (tree, node=NULL, tips=NULL, edges=NULL) {
  # Return all child tips
  .fun <- .phyloOrNodelist (tree, .getParent__phylo,
                            .getParent__NodeList)
  .fun (tree, node=node, tips=tips, edges=edges)
}

#' @name getEdges
#' @title Get edges
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
  .fun <- .phyloOrNodelist (tree, .getEdges__phylo,
                            .getEdges__NodeList)
  .fun (tree, node = node, tips = tips, type = type)
}

#' @name getNodes
#' @title Get parent nodes
#' @description Return all nodes that connect the specified node to the root
#' @details No details
#' @template base_template
#' @param node number indicating node
#' @export
#' @examples
#' # example.var <- exampleFun (test.data)

getNodes <- function (tree, node) {
  .fun <- .phyloOrNodelist (tree, .getNodes__phylo,
                            .getNodes__NodeList)
  .fun (tree, node)
}

#' @name getClades
#' @title Get all children of each node
#' @description All children for each node in the tree are returned.
#' @details Returns a list of lists: 'clade.children' contains all the children of each clade,
#'  'clade.node' contains the clade number of each clade in clade.children.
#' @template base_template
#' @export
#' @examples
#' # example.var <- exampleFun (test.data)

getClades <- function (tree) {
  .fun <- .phyloOrNodelist (tree, .getClades__phylo,
                            .getClades__NodeList)
  .fun (tree)
}

#' @name getNodeLabels
#' @title Get node labels based on online taxonomic database
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
  .fun <- .phyloOrNodelist (tree, .getNodeLabels__phylo,
                            .getNodeLabels__NodeList)
  .fun (tree, all = all, datasource = datasource)
}

#' @name getSubtrees
#' @title Get non-redundant trees of clades from a phylogenetic tree
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
  .fun <- .phyloOrNodelist (tree, .getSubtrees__phylo,
                            .getSubtrees__NodeList)
  .fun (tree, min.n, max.n, verbose = verbose)
}