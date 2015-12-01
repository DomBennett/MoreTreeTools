.addTip__NodeList <- function (tree, edge, tip_name, node_name,
                               node_age, tip_age=0) {
  # edge refers to the sister of the new tip to be added
  sister <- tree@nodelist[[edge]]
  # look up prenode, parent
  parent <- tree@nodelist[[sister$prenode]]
  # key stats
  sister_age <- tree@age - sister$predist
  parent_age <- tree@age - parent$predist
  new_sister_span <- node_age - sister_age
  new_tip_span <- node_age - tip_age
  new_node_pd <- new_tip_span + new_sister_span
  new_node_span <- parent_age - node_age
  new_node_predist <- new_node_span + parent$predist
  new_tip_predist <- new_tip_span + new_node_predist
  # create new tip
  new_tip <- list ('id'=tip_name,
                   'span'=new_tip_span,
                   'prenode'=c (node_name),
                   'postnode'=c (),
                   'children'=c (),
                   'pd'=0,
                   'predist'=new_tip_predist)
  # create new internal node
  new_node <- list ('id'=node_name,
                    'span'=new_node_span,
                    'prenode'=sister$prenode,
                    'postnode'=c (sister$id, tip_name),
                    'children'=c (sister$id, tip_name),
                    'pd'=new_node_pd,
                    'predist'=node_age)
  sister$span <- new_sister_span
  # add to NodeList
  tree@nodelist[[node_name]] <- new_node
  tree@nodelist[[tip_name]] <- new_tip
  # loop throgh each prenode
  node_id <- new_node$prenode
  while (length (node_id) > 0) {
    node <- tree@nodelist[[node_id]]
    node$pd <- node$pd + new_tip_span
    node$children <- c (node$children, tip_name)
    tree@nodelist[[node$id]] <- node
    node_id <- node$prenode
  }
  .update (tree)
}
.removeTip__NodeList <- function (tree, tips, preserve_age) {
  cat ('Not yet implemented for NodeList')
}
.collapseTips__NodeList <- function (tree, min_length, iterative=TRUE) {
  cat ('Not yet implemented for NodeList')
}


