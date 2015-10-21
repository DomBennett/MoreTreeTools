\name{EDBMMtrees}
\alias{EDBMMtrees}
\docType{data}
\title{Simualated phylogenetic trees using EDBMM}
\description{
  List of three large (~1000 tips) simulated \code{phylo} trees with fossils
  and edge diversity data frames generated with \code{calcEdgeDiversity}.
  
  \code{names(EDBMMtrees)} refers to:
  \itemize{
  \item \code{DE} -- Dead-End simulated tree, generated with sig = -1 and eps = +1
  \item \code{Pan} -- Panchronic simulated tree, generated with sig = -1 and eps = -1
  \item \code{Hyd} -- Hydra simulated tree, generated with sig = -1 and eps = 0
  }
  
  Data was generated with this \href{https://github.com/DomBennett/MoreTreeTools/tree/master/other/model_edbmm.R}{script}.
}
\usage{
   data (EDBMMtrees)
}
\format{
  \code{EDBMMtrees} is a \code{list}.
}
\seealso{
  \code{\link{runEDBMM}}, \code{\link{chromatophylo}}
}
\examples{
data (EDBMMtrees) # load object
str (edbmm.trees)  # inspect
# plot each for comparison
chromatophylo (edbmm.trees[['DE']][['tree']],
               edge.cols=edbmm.trees[['DE']][['edge.diversity']],
               legend.title='Diversity') + ggtitle ('DE')
chromatophylo (edbmm.trees[['Pan']][['tree']],
               edge.cols=edbmm.trees[['Pan']][['edge.diversity']],
               legend.title='Diversity') + ggtitle ('Pan')
chromatophylo (edbmm.trees[['Hyd']][['tree']],
               edge.cols=edbmm.trees[['Hyd']][['edge.diversity']],
               legend.title='Diversity') + ggtitle ('Hyd')
}
\keyword{datasets}
