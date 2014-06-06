\name{hominoids}
\alias{hominoids}
\docType{data}
\title{Phylogenetic tree of Hominoids}
\description{
  Phylogenetic tree of the Superfamily Hominoidea, sister group of the Old World Monkeys, trimmed from Bininda-Edmonds et al. (2007)'s Mammalian supertree. Here used for testing and demonstrating MoreTreeTools functions.
}
\usage{
   data (hominoids)
}
\format{
  \code{hominoids} is a \code{phylo} from the pacakge \code{ape}.
}
\source{
  Bininda-Emonds, O. R. P., Cardillo, M., Jones, K. E., MacPhee, R. D. E., Beck, R. M. D., Grenyer, R., … Purvis, A. (2007). The delayed rise of present-day mammals. Nature, 446(7135), 507–12. doi:10.1038/nature05634
}
\examples{
data (hominoids) # load object
plot (hominoids) # plot
}
\keyword{datasets}
