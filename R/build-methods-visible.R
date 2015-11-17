#' @name addTip
#' @title Add a new tip to a phylogenetic tree
#' @description A new tip is added to an existing phylogenetic tree by
#' specifying an edge/branch in the tree, name of the new tip and ages of
#' the new tip and node.
#' @details Return a tree with an extra tip. The new tip is created by
#' finding the \code{edge} in the tree, adding a new node at \code{node.age}
#' and then by adding a new branch to \code{tip.age}.
#' The function is modified code from \code{multi2di} from \code{ape}.
#' N.B. \code{ape} uses node.label to label all internal nodes, \code{MoreTreeTools}
#' uses it to labels all nodes including tip nodes.
#' @template base_template
#' @param edge tree branch where tip will be added
#' @param tip.name name for the new tip to be added
#' @param node.age age of the new node to be created
#' @param tip.age age of the new tip to be added (default 0)
#' @param node.label label for the new
#' @export
#' @examples
#' # example.var <- exampleFun (test.data)

addTip <- function (tree, edge, tip_name, node_name,
                    node_age, tip_age=0) {
  # Safety checks
  if (tip_age > node_age) {
    stop ('node_age must be greater than tip_age')
  }
  if (node_age > getSize (tree, 'rtt')) {
    stop ('node_age given is greater than tree age')
  }
  if (class (tree) == 'phylo') {
    if (tip_name %in% tree$tip.label) {
      stop ('tip_name already in tree')
    }
    if (edge > nrow (tree$edge)) {
      stop ('edge not in tree')
    }
    if (node_name %in% tree$node.label) {
      stop ('node_name already in tree')
    }
  } else if (class (tree) == 'NodeList') {
    if (tip_name %in% tree@tips) {
      stop ('tip_name already in tree')
    }
    if (node_name %in% tree@nodes) {
      stop ('node_name already in tree')
    }
    if (!edge %in% tree@tips) {
      stop ('edge must refer to a tip node ID for a NodeList object')
    }
  } else {
    stop ('tree must be phylo or NodeList')
  }
  .fun <- .phyloOrNodelist (tree, .addTip__phylo,
                            .addTip__NodeList)
  .fun (tree, edge, tip_name, node_name,
        node_age, tip_age)
}

#' @name removeTip
#' @title Remove a tip from a phylogenetic tree
#' @description Remove tips from a phylogenetic tree
#' @details Based on \code{ape}'s \code{drop.tip} but with
#' \code{preserve.age}. This will ensure
#' the root to tip distance of the tree is maintained by adding length
#' to next immediate internal edge whenever the removal of a tip
#' leads to a reduction in a tree's age. If root node is lost, then
#' lost edge length is added to \code{root.edge} slot.
#' @template base_template
#' @param tips vector of tip names
#' @param preserve.age maintain the root to tip distance of the tree
#' @export
#' @examples
#' tree <- rtree (10)
#' par (mfrow=c (1,3))
#' plot (tree)
#' plot (removeTip (tree, tips=c ('t1', 't2', 't3'), preserve.age=FALSE))
#' plot (removeTip (tree, tips=c ('t1', 't2', 't3'), preserve.age=TRUE), root.edge=TRUE)

removeTip <- function (tree, tips, preserve_age) {
  .fun <- .phyloOrNodelist (tree, .removeTip__phylo,
                            .removeTip__NodeList)
  .fun (tree, tips, preserve_age)
}

#' @name collapseTips
#' @title Remove all tips smaller than min.length
#' @description Return a phylo object with tip edges less than min.length
#' removed.
#' @details Tips are removed, but the tree age is preserved. Returned
#' trees are therefore not the same phylogenetic tree as that given.
#' The purpose of the function is to generate representations of
#' trees with fewer tips for plotting purposes.
#' 
#' 
#' Set iterative as TRUE to keep removing all tip edges less than min.length.
#' @template base_template
#' @param min.length smallest length for a tip edge
#' @param iterative run until no tip edges less than min.length, default TRUE
#' @export
#' @examples
#' tree <- rtree (100)
#' par (mfrow = c (1,2))
#' plot (tree, show.tip.label=FALSE, main='Full tree')
#' tree.age <- getSize (tree, 'rtt')
#' ctree <- collapseTips (tree, min.length=tree.age*0.1)
#' plot (ctree, show.tip.label=FALSE, main='Collapsed tree')

collapseTips <- function (tree, min_length, iterative=TRUE) {
  .fun <- .phyloOrNodelist (tree, .removeTip__phylo,
                            .removeTip__NodeList)
  .fun (tree, min_length, iterative)
}