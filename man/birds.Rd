\name{birds}
\Rdversion{1.1}
\alias{birds}
\alias{birdsbreed}
\alias{birdswinter}
\alias{resourceD}
\docType{data}
\title{Birds resource use data}
\description{
  The resource use data for a set of 10 hypothetical bird species during two different seasons. There are five resource classes, and their relationships are given in matrix \code{resourceD}.
}
\usage{data(birds)}
\format{
 birdsbreed - a data frame containing the resource use of 10 birds for a set of five resources during the breeding season.

 birdswinter - a data frame containing the resource use of 10 birds for a set of five resources during the winter season.

 resourceD - a distance matrix (object of class \code{\link{dist}}) containing the relationships between the five resources.
}

\examples{
data(birds)
}
\keyword{datasets}
