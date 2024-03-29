% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pruneindicators.R
\name{pruneindicators}
\alias{pruneindicators}
\title{Determines the best subset of indicators}
\usage{
pruneindicators(
  x,
  At = 0,
  Bt = 0,
  sqrtIVt = 0,
  alpha = 1,
  max.indicators = 4,
  verbose = FALSE
)
}
\arguments{
\item{x}{An object of class '\code{\link{indicators}}'.}

\item{At}{Threshold for positive predictive value. Combinations with lower values are not kept.}

\item{Bt}{Threshold for sensitivity. Combinations with lower values are not kept.}

\item{sqrtIVt}{Threshold for (square root of) indicator value. Combinations with lower values are not kept.}

\item{alpha}{Threshold for statistical significance of indicator value. Combinations with higher p-values are not kept.}

\item{max.indicators}{Maximum number of species combinations to be kept. If \code{NULL}, the function returns all the non-nested valid indicators without further selection.}

\item{verbose}{If TRUE, prints the results of each step.}
}
\value{
An object of class '\code{\link{indicators}}' with only the species combinations selected.
}
\description{
This function allows reducing drastically the number of species combinations to be retained for a given target site group.
}
\details{
First, the function selects those indicators (species or species combinations) with valid positive predictive value, sensitivity and indicator value, according to the input thresholds. If the object '\code{speciescomb}' contains confidence intervals, then the lower bounds are used to select the valid indicators. Second, the function discards those valid indicators whose occurrence pattern is nested within other valid indicators. Third, the function evaluates the \code{\link{coverage}} of the remaining set of indicators and explores subsets of increasing number of indicators, until the same coverage is attained and the set of indicators is returned. If the maximum allowed members is attained (\code{max.indicators}) then the set of indicators with maximum coverage is returned.
}
\examples{
library(stats)

data(wetland) ## Loads species data

## Creates three clusters using kmeans
wetkm <- kmeans(wetland, centers=3) 

## Run indicator analysis with species combinations for the first group
sc <- indicators(X=wetland, cluster=wetkm$cluster, group=1, verbose=TRUE, At=0.5, Bt=0.2)

## Finds the 'best' subset of indicators
sc2 <- pruneindicators(sc, At=0.5, Bt=0.2, verbose=TRUE)
print(sc2)

}
\references{
De \enc{Cáceres}{Caceres}, M., Legendre, P., Wiser, S.K. and Brotons, L. 2012. Using species combinations in indicator analyses. Methods in Ecology and Evolution 3(6): 973-982.

De \enc{Cáceres}{Caceres}, M. and Legendre, P. 2009. Associations between species and groups of sites: indices and statistical inference. Ecology 90(12): 3566-3574.
}
\seealso{
\code{\link{indicators}}, \code{\link{coverage}}
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, EMF-CREAF
}
