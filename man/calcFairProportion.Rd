\name{calcFairProportion}
\alias{calcFairProportion}
\title{Calculate Evolutionary Distinctiveness (Fair Proportion)}
\usage{
calcFairProportion(tree, tips = "all", progress = "none")
}
\arguments{
  \item{tree}{phylogenetic tree (ape class)}

  \item{tips}{vector of tips for which FP is calculated,
  else 'all'}

  \item{progress}{bar showing the progress of the
  calculation. Either time, text, tk or none. Uses
  \code{plyr}'s \code{create_progress_bar}.}
}
\description{
Calculate the proportion of a phylogenetic tree that
species represent
}
\details{
Calculate the evolutionary distinctivess of species using
Issac et al. (2007)'s Fair Proportion.
}
\examples{
# example.var <- exampleFun (test.data)
}
\references{
No references yet!
\url{https://github.com/DomBennett/MoreTreeTools}.
}
\keyword{ecology}
\keyword{evolution}
\keyword{manipulation}
\keyword{palaeontology}
\keyword{phylogeny}

