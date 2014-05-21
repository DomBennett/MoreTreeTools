\name{reduceTree}
\alias{reduceTree}
\title{Reduce tree by taxonomic rank}
\usage{
reduceTree(tree, level, datasource)
}
\arguments{
  \item{tree}{phylogenetic tree (ape class)}

  \item{level}{rank by which tips will be reduce (e.g.
  kingdom)}

  \item{datasource}{a number indicating the GNR datasource
  to search against}
}
\description{
Reduce tips by identifying tips in the same taxonomic group
at a specifed taxonomic rank
}
\details{
Search Global Names Resolver to find shared taxonomic
groups of all tips. Reduce tree by condensing members of
the same taxonomic group into a single tip. The user can
specify at what rank the tree should be reduced. N.B.
unresolved names are dropped from the tree and are not
represented in tip counts.
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
