# Thinking about creating a class for faster tree manipulation

# TODO:
# -- create addNode
# -- create removeNode
# -- create updateNodeList (updates values given a change)
# -- create converter from Phylo to NodeList

# Process for creating a class
# 1. Use setClass to define the class
# 2. Use setGeneric to define a method for the class
# 3. Use setMethod to define th function for the method

# Node declaration
node_representation <- representation (children='vector',
                                       min_age='numeric',
                                       max_age='numeric',
                                       prev_node='integer',
                                       next_node='integer',
                                       node_id='integer')
setClass ('Node', representation=node_representation)
setGeneric ('print')
setMethod ('print', c('x'='Node'),
           function(x){
             cat ('Node Object:\n\tID = [', x@node_id, '], children = [',
                  length (x@children), '], age = [',
                  x@max_age, ']', sep='')
           })
setMethod ('show', 'Node',
           function(object){
             print (object)
           })
# NodeList declaration
nodelist_representation <- representation (nodes='list')
prototype_root_node <- new ('Node', node_id=1L)
nodelist_prototype <- prototype (nodes=list(prototype_root_node))
setClass ('NodeList', representation=nodelist_representation,
          prototype=nodelist_prototype)
setGeneric ('print')
setMethod ('print', c('x'='NodeList'),
           function(x){
             cat ('NodeList Object of [',
                  length (x@nodes), '] nodes', sep='')
           })
setMethod ('show', 'NodeList',
           function(object){
             print (object)
           })
setMethod ('[', c ('NodeList', 'numeric', 'missing', 'ANY'),
           function(x, i, j, ..., drop=TRUE) {
             i <- as.integer (i)
             initialize(x, nodes=x@nodes[i])
           })
setMethod ('[[', c ('NodeList', 'numeric', 'missing'),
           function(x, i, j, ...) {
             i <- as.integer (i)
             x@nodes[[i]]
           })
setGeneric ('addNode', signature= c('x', 'parent_node', 'min_age',
                                    'max_age', 'tip'),
            function(x, parent_node, min_age, max_age, tip) {
              genericFunction ('addNode')
              })
setMethod ('addNode', c('x'='NodeList', 'parent_node'='integer',
                        'min_age'='numeric', 'max_age'='numeric',
                        'tip'=ANY),
           function(x, parent_node, min_age, max_age, tip=NA) {
             node_id <- as.integer (length (x@nodes) + 1)
             children <- c (x@nodes[[parent_node]]@children, tip)
             node <- new ('Node', min_age=min_age, max_age=max_age,
                          prev_node=parent_node, node_id=node_id,
                          children=children)
             x@nodes[[node_id]] <- node
             return (x)
           })
nodelist <- new('NodeList')
nodelist <- addNode(nodelist, parent_node=1L, min_age=0, max_age=1, tip=NA)
nodelist[3:4]
str(nodelist@nodes[[2]])

initialize (nodelist, nodes=list (nodelist@nodes[[1]]))

# Tree declaration
tree_representation <- representation (nodes='NodeList')
setClass ('TreeList', representation=tree_representation,
          prototype=prototype (nodes=new('NodeList')))
new ('TreeList')