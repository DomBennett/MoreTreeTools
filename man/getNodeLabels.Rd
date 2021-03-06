% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get-methods.R
\name{getNodeLabels}
\alias{getNodeLabels}
\title{Get node labels based on online taxonomic database}
\usage{
getNodeLabels(tree, all = FALSE, cache = FALSE, parent = NULL)
}
\arguments{
\item{tree}{phylogenetic tree (ape class)}

\item{all}{count tips as nodes, default False}

\item{cache}{T/F, create a local cache of downloaded names?}

\item{parent}{specify parent of all names to prevent false names}
}
\description{
Return names of each node in tree based on searching tip labels
 through Global Names Resolver \url{http://resolver.globalnames.org/} in NCBI.
}
\details{
For each node, all the children are searched, the taxonomic lineages returned and
then searched to find the lowest shared name.
All the tip labels are searched against a specified taxonomic database through the GNR and NCBI.
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

