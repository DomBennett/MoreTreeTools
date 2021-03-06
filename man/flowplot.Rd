% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/experimental-methods.R
\name{flowplot}
\alias{flowplot}
\title{Plot tree as turbulence flow image (experimental)}
\usage{
flowplot(tree, n = getSize(tree))
}
\arguments{
\item{tree}{phylogenetic tree (ape class)}

\item{n}{number of pixels}
}
\description{
Experimental plot where a tree is broken down into
X and Y pixels. Edges are mapped to this Cartesian space with
colour representing number of descendents. The idea is to better
visualise very large trees, particular ones with fossils. Currently only
works with ultrametric trees.
}
\details{
No details
}
\examples{
# Generate and plot a random tree
tree <- compute.brlen (rtree (100))
flowplot (tree)
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

