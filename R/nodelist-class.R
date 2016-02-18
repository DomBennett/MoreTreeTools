#' @name NodeList
#' @title NodeList Class
#' @description \code{MoreTreeTool}'s S4 Class for representing phylogenetic trees.
#' \code{NodeList} objects hold a list of \code{Node} Reference Class Objects. Each
#' \code{Node} points to connecting nodes.
#' @exportClass NodeList
#' @examples
#' # No examples yet
# TODO: create validity check
setClass ('NodeList', representation=representation (
  nodelist='list',       # list of Node objects
  nodes='vector',        # vector of Node ids that are internal nodes
  tips='vector',         # vector of Node ids that are tips
  age='numeric',         # numeric of max root to tip distance
  pd='numeric',          # numeric of total branch length of tree
  extant='vector',       # vector of Node ids of all tips with 0 age
  extinct='vector',      # vector of Node ids of all tips with age > 0
  brnchlngth='logical',  # logical, do nodes have span
  ultrmtrc='logical',    # logical, do all tips end at 0
  plytms='logical',      # logical, is tree bifurcating
  tol='numeric',         # numeric of tolerance for determining extant
  root='character'),     # character of Node id of root
  prototype=prototype (tol=1e-8))

# Display methods
setMethod ('as.character', c('x'='NodeList'),
           function(x) {
             paste0 ('Tree (NodeList Object) of [', length (x@tips),'] tips')
           })
setMethod ('show', 'NodeList',
           function(object){
             print (object)
           })
setMethod ('str', c('object'='NodeList'),
           function (object, max.level=2L, ...) {
             str@default (object, max.level=max.level, ...)
           })
setGeneric ('print')
setMethod ('print', c('x'='NodeList'),
           function(x){
             msg <- 'Tree (NodeList Object):\n'
             msg <- paste0 (msg, '  -- [', length (x@tips), '] tips\n')
             msg <- paste0 (msg, '  -- [', length (x@nodes), '] internal nodes\n')
             if (x@plytms) {
               msg <- paste0 (msg, '  -- Polytomous\n')
             } else {
               msg <- paste0 (msg, '  -- Binary\n')
             }
             if (is.na (rootNode (x))) {
               if (is.na (age (x))) {
                 msg <- paste0 (msg, '  -- Unrooted and without branch lengths\n')
               } else {
                 msg <- paste0 (msg, '  -- Unrooted, with branch lengths\n')
                 msg <- paste0 (msg, '  -- PD [', signif (x@pd, 3), ']\n')
               }
             } else {
               msg <- paste0 (msg, '  -- Age [', signif (x@age, 3), ']\n')
               msg <- paste0 (msg, '  -- PD [', signif (x@pd, 3), ']\n')
               if (x@ultrmtrc) {
                 msg <- paste0 (msg, '  -- Ultrametric (all tips are extant)\n')
               } else {
                 msg <- paste0 (msg, '  -- Not ultrametric (with extinct tips)\n')
               }
               msg <- paste0 (msg, '  -- Root node is ["', x@root, '"]\n')
             }
             cat (msg)
           })

# Manip methods
setMethod ('[[', c ('NodeList', 'character', 'missing'),
           function(x, i, j, ...) {
             x@nodelist[[i]]
           })
setGeneric ('add', signature=c('x', 'node'),
            function(x, node) {
              genericFunction ('add')
            })
setMethod ('add', c ('NodeList', 'Node'),
           function (x, node) {
             x@nodelist[[node$id]] <- node
             .update(x)
           })
setGeneric ('remove', signature=c('x', 'node'),
            function(x, node) {
              genericFunction ('remove')
            })
setMethod ('remove', c ('NodeList', 'Node'),
           function (x, node) {
             for (n in x@nodelist[[node$id]]$postnode) {
               remove (x, n)
             }
             # remove from prenodes
             for (n in node$prenode) {
               pids <- sapply (n$postnode, function (x) x$id)
               n$postnode <- n$postnode[pids != node$id]
             }
             # remove from nodelist
             x@nodelist <- x@nodelist[names (x@nodelist) != node$id]
             #gc (verbose=FALSE)
             .update(x)
           })
setGeneric ('.update', signature=c('x'),
            function(x) {
              genericFunction ('.update')
            })
setMethod ('.update', 'NodeList',
            function (x) {
              with_pstndes <- sapply (x@nodelist,
                                      function (x) length (x$postnode) == 0)
              x@tips <- names (with_pstndes)[with_pstndes]
              x@nodes <- names (with_pstndes)[!with_pstndes]
              x@brnchlngth <- all (sapply (x@nodelist, function (x) length (x$span) > 0))
              if (x@brnchlngth) {
                if (length (x@root) > 0) {
                  x@age <- x@nodelist[[x@root]]$postDist()
                  extant_is <- unlist (sapply (x@tips, function (i) {
                    length (x@nodelist[[i]]$postnode) == 0 &
                      (x@age - x@nodelist[[i]]$preDist()) <= x@tol}))
                  x@extant <- names (extant_is)[extant_is]
                  x@extinct <- x@tips[!x@tips %in% x@extant]
                  x@ultrmtrc <- all (x@tips %in% extant (x))
                }
                x@pd <- sum (unlist (sapply (x@nodelist, function (x) x$span)))
              }
              x@plytms <- any (sapply (x@nodelist, function (x) length (x$postnode) > 2))
              initialize (x)
            })

# Accessor method
setGeneric ('setTol', signature=c('x', 'n'),
            function(x, n) {
              genericFunction ('setTol')
            })
setMethod ('setTol', c ('NodeList', 'numeric'),
           function(x, n){
             x@tol <- n
             .update (x)
           })

# Info methods (user friendly)
setGeneric ('nTips', signature=c('x'),
            function(x) {
              genericFunction ('nTips')
            })
setMethod ('nTips', 'NodeList',
           function(x){
             length (x@tips)
           })
setGeneric ('plytms', signature=c('x'),
            function(x) {
              genericFunction ('plytms')
            }) 
setMethod ('plytms', 'NodeList',
           function(x) {
             x@plytms
           })
setGeneric ('ultrmtrc', signature=c('x'),
            function(x) {
              genericFunction ('ultrmtrc')
            }) 
setMethod ('ultrmtrc', 'NodeList',
           function(x) {
             x@ultrmtrc
           })
setGeneric ('extant', signature=c('x'),
            function(x) {
              genericFunction ('extant')
            }) 
setMethod ('extant', 'NodeList',
           function(x) {
             if (length (x@extant) == 0) {
               return (NULL)
             }
             x@extant
           })
setGeneric ('extinct', signature=c('x'),
            function(x) {
              genericFunction ('extinct')
            }) 
setMethod ('extinct', 'NodeList',
           function(x) {
             if (length (x@extinct) == 0) {
               return (NULL)
             }
             x@extinct
           })
setGeneric ('nNodes', signature=c('x'),
            function(x) {
              genericFunction ('nNodes')
            })
setMethod ('nNodes', 'NodeList',
           function(x){
             length (x@nodes)
           })
setGeneric ('rootNode', signature=c('x'),
            function(x) {
              genericFunction ('rootNode')
            })
setMethod ('rootNode', 'NodeList',
           function(x) {
             if (length (x@root) == 0) {
               return (NA)
             }
             x@root
           })
setGeneric ('age', signature=c('x'),
            function(x) {
              genericFunction ('age')
            })
setMethod ('age', 'NodeList',
           function(x){
             if (length (x@age) == 0) {
               return (NA)
             }
             x@age
           })
setGeneric ('pd', signature=c('x'),
            function(x) {
              genericFunction ('pd')
            })
setMethod ('pd', 'NodeList',
           function(x){
             if (length (x@pd) == 0) {
               return (NA)
             }
             x@pd
           })