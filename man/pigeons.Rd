\name{pigeons}
\Rdversion{1.1}
\alias{pigeons}
\alias{dfood}
\alias{diet.barcelona}
\alias{diet.moia}
\docType{data}
\title{Pigeon seed consumption data}
\description{
  Data set used for the dietary analysis of two feral pigeon (Columba livia) populations. There are six seed types, and their distance relationships are given in matrix \code{dfood}.
}
\usage{data(pigeons)}
\format{
 diet.barcelona - a data frame containing the seed consumption (in grams) of pigeons taken from the population in Barcelona.

 diet.moia - a data frame containing the seed consumption (in grams) of pigeons taken from the population in Moia.

 dfood - a distance matrix (object of class \code{\link{dist}}) containing the relationships between the six seed types. Distances were calculated from data describing the nutritional content and dimensions of each seed type.
}
\references{
De Caceres, M., Sol, D., Lapiedra, O. and P. Legendre. (2011) A framework for estimating niche metrics using the resemblance between qualitative resources. Oikos 120: 1341-1350.
}
\examples{
data(pigeons)
}
\keyword{datasets}
