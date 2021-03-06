% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get-methods.R
\name{getEdges}
\alias{getEdges}
\title{Get edges}
\usage{
getEdges(tree, node = NULL, tips = NULL, type = 1)
}
\arguments{
\item{tree}{phylogenetic tree (ape class)}

\item{node}{number indicating node}

\item{tips}{tip labels or numbers}

\item{type}{1, 2 or 3 indicating which method to use}
}
\description{
Return all children edges from a node,
 or all edges connected to tips based on three different methods
}
\details{
If node specified, all edges that descend from specified node are returned

Otherwise, if tips are specified edges are returned based on one of three methods:
Type:

1 -- phylogeny consisting solely of the tips, default

2 -- edges from taxon tips to terminal node

3 -- edges unique to tips
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

