.addTip__NodeList <- function (tree, edge, tip_name, node_name,
                               node_age, tip_age=0) {
  # edge refers to the sister tip
  sister <- tree@nodelist[[edge]]
  # look up prenode and up date
  parent <- sister$prenode[[1]]
  parent_age <- tree@age - parent$preDist()
  # create new internal node
  new_node <- new ('Node', id=node_name, span=parent_age-node_age,
                   prenode=c (parent))
  # create new tip
  new_tip <- new ('Node', id=tip_name, span=node_age-tip_age,
                  prenode=c (new_node))
  new_node$postnode <- c (sister, new_tip)
  # update parent's postnode
  parent_postnodes <- parent$postnode[which (sapply (parent$postnode,
                                                     function (x) x$id) != edge)]
  parent_postnodes <- c (parent_postnodes, new_node)
  # update sister's prenode
  sister$prenode <- c (new_node)
  # add to NodeList
  tree@nodelist[[node_name]] <- new_node
  tree@nodelist[[tip_name]] <- new_tip
  # Update
  .update (tree)
}
.removeTip__NodeList <- function (tree, tips, preserve_age) {
  cat ('Not yet implemented for NodeList')
}
.collapseTips__NodeList <- function (tree, min_length, iterative=TRUE) {
  cat ('Not yet implemented for NodeList')
}


