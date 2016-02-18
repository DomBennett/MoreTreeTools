.addTip__NodeList <- function (tree, edge, tip_name, node_name,
                               node_age, tip_age=0) {
  # edge refers to the sister of the new tip to be added
  sister <- tree@nodelist[[edge]]
  # look up prenode and up date
  parent <- sister$prenode[[1]]
  parent_age <- .getAge__NodeList (tree, parent)
  # create new internal node
  new_node <- list ('id'=node_name, 'span'=parent_age-node_age,
                   'prenode'=c (parent))
  # create new tip
  new_tip <- list ('id'=tip_name, 'span'=node_age-tip_age,
                  'prenode'=c (new_node$id))
  new_node$postnode <- c (sister$id, new_tip$id)
  # update parent's postnode
  parent_postnodes <- tree@nodelist[[parent]]$postnode
  parent_postnodes <- parent_postnodes[parent_postnodes != edge]
  tree@nodelist[[parent]]$postnode <- c (parent_postnodes, new_node$id)
  # update sister's prenode
  tree@nodelist[[edge]]$prenode <- c (new_node$id)
  # add to NodeList
  tree@nodelist[[node_name]] <- new_node
  tree@nodelist[[tip_name]] <- new_tip
  # Update
  tree
}
.removeTip__NodeList <- function (tree, tips, preserve_age) {
  cat ('Not yet implemented for NodeList')
}
.collapseTips__NodeList <- function (tree, min_length, iterative=TRUE) {
  cat ('Not yet implemented for NodeList')
}


