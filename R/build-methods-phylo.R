.addTip__phylo <- function (tree, edge, tip_name, node_name,
                            node_age, tip_age=0) {
  insert <- function (target_vector, source_vector, index) {
    #http://stackoverflow.com/questions/1493969/how-to-insert-elements-into-a-vector
    index <- index - 1
    positions <- c (seq_along (target_vector), index + 0.5)
    res <- c (target_vector, source_vector)
    res[order (positions)]
  }
  # nodes that make edge
  node_1 <- tree$edge[edge, 1]
  node_2 <- tree$edge[edge, 2]
  new_tip_edge_length <- node_age - tip_age
  if (new_tip_edge_length < 0) {
    stop ("Node age must be greater than tip age.")
  }
  if (!is.null (tree$node_ages)) {
    tree_node_ages <- tree$node_ages
    node_1_age <- tree_node_ages[node_1]
  } else {
    node_1_age <- getAge (tree, node=node_1)$age
  }
  new_node_edge_length <- node_1_age - node_age
  if (new_node_edge_length < 0) {
    stop ("Node age is greater than edge.")
  }
  # build a new edge matrix
  edges_to_replace <- c(edge, getEdges (tree, node_2))
  new_edge_lengths <- c (new_node_edge_length, new_tip_edge_length,
                         tree$edge.length[edge] - new_node_edge_length,
                         tree$edge.length[edges_to_replace][-1])
  target <- node_1 + 1
  tree$edge[tree$edge > length (tree$tip.label)] <- 
    tree$edge[tree$edge > length (tree$tip.label)] + 1
  tree$Nnode <- tree$Nnode + 1
  new_tip_labels <- c (tree$tip.label, tip_name)
  new_tip <- length (new_tip_labels)
  new_node <- new_tip + tree$Nnode
  new_edges <- matrix (nrow = length (edges_to_replace) + 2, ncol = 2)
  new_edges[1, ] <- c (target, new_node)
  new_edges[2, ] <- c (new_node, new_tip)
  for (i in 1:length (edges_to_replace)) {
    new_edge <- tree$edge[edges_to_replace[i],]
    new_edge[new_edge == target] <- new_node
    new_edges[i+2, ] <- new_edge
  }
  # replace old with new
  tree$edge <- rbind (tree$edge[-edges_to_replace, ], new_edges)
  tree$edge.length <- c (tree$edge.length[-edges_to_replace],
                         new_edge_lengths)
  tree$tip.label <- new_tip_labels
  # tidy up new tree (multi2di code)
  if (!is.null (attr (tree, "order"))) 
    attr(tree, "order") <- NULL
  if (!is.null (tree$node.label)) {
    tree$node.label <- insert (tree$node.label, node_name,
                               new_node - length (tree$tip.label))
  }
  if (!is.null (tree$node_ages)) {
    tree_node_ages <- insert (tree_node_ages, c (tip_age, node_age),
                              c (new_tip, new_node))
  }
  tree <- reorder (tree)
  newNb <- integer (tree$Nnode)
  n <- length (tree$tip.label)
  newNb[1] <- n + 1L
  sndcol <- tree$edge[, 2] > n
  o <- 1 + rank (tree$edge[sndcol, 2])
  int_node_i <- (length (tree$tip.label) + 1):
    (length (tree$tip.label) + tree$Nnode)
  if (!is.null (tree$node.label)) {
    o <- 1 + rank (tree$edge[sndcol, 2])
    tree$node.label <- tree$node.label[c (1, o)]
  }
  if (!is.null (tree$node_ages)) {
    int_node_age <- tree_node_ages[int_node_i]
    int_node_age <- int_node_age[c(1, o)]
    tree$node_ages[int_node_i] <- int_node_age
  }
  tree$edge[sndcol, 2] <- newNb[tree$edge[sndcol, 2] - n] <-
    n + 2:tree$Nnode
  tree$edge[, 1] <- newNb[tree$edge[, 1] - n]
  tree
}
.removeTip__phylo <- function (tree, tips, preserve_age) {
  # internal
  .run <- function (tip_name) {
    # names description:
    # - tree_age, age of the tree from root to tip
    # - root node, the node number of the root of the tree
    # - max_node, the maximum internal node number of the tree
    # - edges_to_drop, edge numbers that will be dropped (max 2, target and parent edges)
    # - target_node, node number of tip_name
    # - target_edge, connecting edge of tip node
    # - node_1, first internal node connecting to target_node and sister_node
    # - node_2, second internal node connecting to tree to node_1
    # - sister_node, node of sister
    # - sister_edge, edge of sister
    # - parent_edge, edge of node_2 and node_1
    # - _ori, original data provided
    # - _mod, modified without missing elements
    # - _new, modified with missing elements 
    #print (tip_name)
    # base vars
    tree_age <- getSize (tree, 'rtt')
    edges_ori <- edges_mod <- tree$edge
    root_node <- length (tree$tip.label) + 1
    max_node <- tree$Nnode + length (tree$tip.label)
    target_node <- which (tree$tip.label == tip_name)
    target_edge <- which (tree$edge[ ,2] == target_node)
    edges_to_drop <- target_edge
    node_1 <- tree$edge[target_edge, 1]
    # get sister
    sister_edge <- which (tree$edge[, 1] == node_1)
    sister_edge <- sister_edge[sister_edge != target_edge]
    sister_node <- tree$edge[sister_edge, 2]
    if (node_1 != root_node) {
      # get parent
      parent_edge <- which (tree$edge[ ,2] == node_1)
      node_2 <- tree$edge[parent_edge,1]
      edges_to_drop <- c (edges_to_drop, parent_edge)
      # update length of sister_edge
      parent_edge_length <- tree$edge.length[parent_edge]
      sister_edge_length <- tree$edge.length[sister_edge]
      target_edge_length <- tree$edge.length[target_edge]
      sister_edge_length_new <- sister_edge_length +
        parent_edge_length
      if (target_edge_length > sister_edge_length &
            preserve_age) {
        sister_edge_length_new <- sister_edge_length_new +
          (target_edge_length - sister_edge_length)
      }
      tree$edge.length[sister_edge] <- sister_edge_length_new
      # replace node_1 with node_2 for sister_edge
      edges_mod[sister_edge, 1] <- node_2
      # remove edges to drop
      edges_new <- edges_mod[-edges_to_drop, ]
      # re-number nodes
      # remove 1 from all nodes greater than or equal to
      #  the node_1
      edges_new[edges_new[ ,1] >= node_1, 1] <-
        edges_new[edges_new[ ,1] >= node_1, 1] - 1
      edges_new[edges_new[ ,2] >= node_1, 2] <-
        edges_new[edges_new[ ,2] >= node_1, 2] - 1
      # remove 1 from all nodes above tip node
      edges_new[edges_new[ ,1] > target_node,1] <-
        edges_new[edges_new[ ,1] > target_node,1] - 1
      edges_new[edges_new[ ,2] > target_node,2] <-
        edges_new[edges_new[ ,2] > target_node,2] - 1
    } else {
      # if node is root, there is no parent edge, must remove sister
      edges_to_drop <- c (edges_to_drop, sister_edge)
      # remove edges to drop
      edges_new <- tree$edge[-edges_to_drop, ]
      # re-number nodes
      edges_new <- edges_new - 1
      edges_new[ ,1] <- edges_new[ ,1] - 1
      edges_new[edges_new[ ,2] >= node_1, 2] <-
        edges_new[edges_new[ ,2] >= node_1, 2] - 1
      if (preserve_age) {
        sister_node_age <- getAge (tree, node=sister_node)[ ,'age']
        tree$root.edge <- tree_age - sister_node_age
      }
    }
    # replace old with new
    tree$edge <- edges_new
    tree$edge.length <- tree$edge.length[-edges_to_drop]
    if (!is.null (tree$edge.label)) {
      tree$edge.label <- tree$edge.label[-edges_to_drop]  
    }
    tree$tip.label <- tree$tip.label[-target_node]
    # update Nnode
    tree$Nnode <- tree$Nnode - 1
    # tidy up new tree
    if (!is.null (attr (tree, "order"))) {
      attr(tree, "order") <- NULL
    }
    #storage.mode (tree$edge) <- "integer"
    #tree_c <- collapse.singles (tree)
    tree <<- tree
  }
  if (!all (tips %in% tree$tip.label)) {
    stop ('Not all tips are in the tree')
  }
  if (length (tips) >= (length (tree$tip.label) - 2)) {
    stop ('Removing too many tips, smallest resulting tree is 3 tips')
  }
  loop_data <- data.frame (tip_name=tips, stringsAsFactors=FALSE)
  m_ply (.data=loop_data, .fun=.run)
  tree
}
.collapseTips__phylo <- function (tree, min_length, iterative=TRUE) {
  while (TRUE) {
    tip_edges <- tree$edge[ ,2] < getSize (tree)
    drop_tip_edges <- (tree$edge.length < min_length) &
      tip_edges
    to_drop <- tree$tip.label[tree$edge[drop_tip_edges,2]]
    if (length (to_drop) == 0) {
      break
    }
    tree <- try (removeTip (tree, tip=to_drop, preserve_age=TRUE),
                 silent=TRUE)
    if (grepl ('Removing too many tips', tree[[1]][1])) {
      stop ('min_length is too great, too many tips are being dropped')
    }
    if (!iterative) {
      break
    }
  }
  tree
}