\name{EDBMMtrees}
\alias{EDBMMtrees}
\docType{data}
\title{Simualated phylogenetic trees using EDBMM}
\description{
  List of three large (~1000 tips) simulated \code{phylo} trees with fossils
  and edge diversity data frames generated with \code{calcEdgeDiversity}.
  
  \code{names(EDBMMtrees)} refers to:
  
  DE -- Dead-End simulated tree, generated with sig = -1 and eps = +1
  Pan -- Panchronic simulated tree, generated with sig = -1 and eps = -1
  Hyd -- Hydra simulated tree, generated with sig = -1 and eps = 0
}
\usage{
   data (EDBMMtrees)
}
\format{
  \code{EDBMMtrees} is a \code{list}.
}
\examples{
data (EDBMMtrees) # load object
str (edbmm.trees)  # inspect
# plot each for comparison
treeplot (edbmm.trees[['DE']][['tree']],
          edge.cols=edbmm.trees[['DE']][['edge.diversity']],
          legend.title='Diversity') + ggtitle ('DE')
treeplot (edbmm.trees[['Pan']][['tree']],
          edge.cols=edbmm.trees[['Pan']][['edge.diversity']],
          legend.title='Diversity') + ggtitle ('Pan')
treeplot (edbmm.trees[['Hyd']][['tree']],
          edge.cols=edbmm.trees[['Hyd']][['edge.diversity']],
          legend.title='Diversity') + ggtitle ('Hyd')
}
\keyword{datasets}
