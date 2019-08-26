\encoding{UTF-8}
\name{combinespecies}
\alias{combinespecies}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Combines species from a community table}
\description{
  Creates species combinations to be used in indicator value analyses. 
}
\usage{
combinespecies(X, min.order = 1, max.order = 3, min.occ = 1, FUN = min, verbose = FALSE, 
add.names = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{A community data table with sites in rows and species in columns. This table can contain either presence-absence or abundance data.}
  \item{min.order, max.order}{Minimum and maximum number of species conforming species combinations.}
  \item{min.occ}{Threshold for minimum occurrence. Combinations with lower values are not kept.}
  \item{FUN}{Function to be calculated for all species forming the species combination. Accepted values are \code{min}, \code{max}, \code{sum} and \code{mean}.}
  \item{verbose}{If TRUE, prints the results of each step.}
  \item{add.names}{If TRUE, adds the names of the species combinations to the output matrix. Species combination names are lists of species concatenated using character '+'.}
  \item{...}{Additional arguments for function \code{FUN}.}
}
\details{
   This function allows creating a data table where rows are sites and columns are combinations of species. Values for a given column of this matrix are derived from the abundance values of the species forming the corresponding combination. In particular, the abundance value for a given combination in a given site is equal to the value returned by function 'FUN' (normally the minimum) among the site values of all species forming the combination. The matrix 'XC' returned by this function can be used in functions \code{\link{strassoc}} and \code{\link{signassoc}}. Alternatively, \code{\link{indicators}} and related functions provide a more elaborated way to explore the indicator value of the simultaneous occurrence of sets of species (i.e. species combinations). 
}
\value{
  An list with:
  \item{XC}{A matrix containing the abundance/occurrence of each species combination.}
  \item{C}{A binary matrix describing the set of species forming each combination ('0' means that the species is not included, and '1' means that the species is included).}
}
\references{
De \enc{Cáceres}{Caceres}, M., Legendre, P., Wiser, S.K. & Brotons, L (in preparation). Using species combinations in indicator analyses.
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CTFC
}
\seealso{\code{\link{indicators}}, \code{\link{strassoc}}, \code{\link{signassoc}}}
\examples{
library(stats)

data(wetland) ## Loads species data

## Create species combinations
Y=combinespecies(X=wetland, max.order=3, min.occ=5, verbose=TRUE)

## Creates three site groups using kmeans
wetkm = kmeans(wetland, centers=3) 

## Calculate indicator value of species combinations for each of the three site groups
strassoc(Y$XC, cluster=wetkm$cluster,func="IndVal.g") 

## Calculate point biserial correlation value of species combinations 
## for each of the three site groups
strassoc(Y$XC, cluster=wetkm$cluster,func="r.g") 

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.