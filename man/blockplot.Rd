\name{blockplot}
\alias{blockplot}
\title{Plot Blocks for Detecting Phylogenetic Signal}
\usage{
blockplot(tree, trait, title = NULL)
}
\arguments{
  \item{tree}{phylogenetic tree (ape class)}

  \item{trait}{vector or matrix of traits for each tip}

  \item{title}{title of plot (default NULL)}
}
\description{
A plot function for visualing phylogenetic signal in a tree
using blocks to represent branches.
}
\details{
No details
}
\examples{
# Generate a balanced tree
tree <- compute.brlen (stree (16, 'balanced'))
# Generate perfect phylogenetic signal traits
# (Note: trait vector should be named)
trait <- c (rep (1, 8), rep (0, 8))
names (trait) <- tree$tip.label
# Plot blocks ...
blockplot (tree, trait, 'Perfect signal')
# Repeat with evenly distributed traits ...
trait <- rep (c (0, 1), 8)
names (trait) <- tree$tip.label
blockplot (tree, trait, 'Maximally overdispered')
# The more evenly distributed, the more symetrical
# the pattern
# A random distribution of traits ...
trait <- randCommData (tree, nsites = 1, nspp = 8)[1, ]
blockplot (tree, trait, 'Random trait')
# ... isn't as symetrical
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

