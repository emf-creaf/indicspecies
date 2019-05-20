\encoding{UTF-8}
\name{multipatt}
\alias{multipatt}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Multi-level pattern analysis}
\description{
 This function studies the association between species patterns and combinations of groups of sites. 
}
\usage{
multipatt(x, cluster, func = "IndVal.g", duleg=FALSE, restcomb=NULL, 
         min.order = 1, max.order = NULL, control=how(), 
         print.perm=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Community data table}
  \item{cluster}{A vector representing a partition of sites}
  \item{func}{Species-site group association function. Four values are accepted \code{"IndVal"}, \code{"IndVal.g"}, \code{"r"} and \code{"r.g"} (lowercase values are also accepted).}
  \item{duleg}{If TRUE, site group combinations are not considered, only the original site groups, like in Dufrêne & Legendre (1997). Internally, \code{duleg = TRUE} equals \code{max.order = 1}.}
  \item{restcomb}{A vector of integer values used to restrict the combinations of site groups to those with ecological sense according to the analyst. The default \code{NULL} indicates that all combinations are used. If \code{duleg=TRUE} this argument is ignored.} 
  \item{min.order}{An integer indicating the minimum order of site group combinations (by default \code{max.order=1} for singletons). Cannot be larger than \code{max.order}.}
  \item{max.order}{An integer indicating the maximum order of site group combinations to be considered: \code{max.order=1} for singletons, \code{max.order=2} for pairs, \code{max.order=3} for triplets... As \code{restcomb}, this parameter provide a way to restrict the site group combinations that make ecological sense. By default all possible site group combinations are considered. If \code{max.order=1} then the function will behave as if \code{duleg=TRUE}.}
  \item{control}{a list of control values describing properties of the permutation design, as returned by a call to \code{\link{how}}.}
  \item{print.perm}{If TRUE, prints permutation numbers after each set of 100 permutations.}
}
\details{
   This function creates combinations of the input clusters and compares each combination with the species in the input matrix x. For each species it chooses the combination with a highest association value. Best matching patterns are tested for statistical significance of the associations. Four association indices are possible (some less than for \code{\link{strassoc}}): "IndVal", "IndVal.g", "r" and "r.g". Indicator value indices will return the pattern that better matches the species observed pattern, whereas correlation indices will return the pattern that creates a highest inside/outside difference. Details are given in De \enc{Cáceres}{Caceres} et al. (2010). The user can restrict the combinations in three ways: (1) by using \code{duleg=TRUE}, which leads to consider single site-groups only; (2) by setting the minimum and maximum order of combinations using \code{min.order} and \code{max.order}; or (3) by using \code{restcomb} to restrict combinations at will. In order to carry out the third way, values in \code{restcomb} must be the indices of combinations that appear in the column \code{index} of object \code{sign} (see below). Complex permutation designs are allowed through the function \code{\link{how}} from package "permute".
}
\value{
  An object of class \code{multipatt} with:
  \item{func}{The name of the function used.}
  \item{comb}{A matrix describing all the combinations studied.}
  \item{str}{A matrix the association strength for all combinations studied.}
  \item{A}{If \code{func = "IndVal"} (or \code{func = "IndVal.g"}) a matrix whose values are the "A" (or "A.g") component of indicator values. Otherwise this element is left as \code{NULL}.}
  \item{B}{If \code{func = "IndVal"} (or \code{func = "IndVal.g"}) a matrix whose values are the "B" component of indicator values. Otherwise this element is left as \code{NULL}.}
  \item{sign}{Data table with results of the best matching pattern, the association value and the degree of statistical significance of the association (i.e. p-values from permutation test). Note that p-values are not corrected for multiple testing.}
}
\note{
This function gives the same results as function \code{indval} in package "labdsv" when used setting \code{func="IndVal.g"} and \code{duleg=TRUE}, excepting the fact that the square root IndVal values is returned instead of the original IndVal.
}
\references{
De \enc{Cáceres}{Caceres}, M. and Legendre, P. 2009. Associations between species and groups of sites: indices and statistical inference. Ecology 90(12): 3566-3574.

De \enc{Cáceres}{Caceres}, M., Legendre, P., Moretti, M. 2010. Improving indicator species analysis by combining groups of sites. Oikos 119(10): 1674-1684.

\enc{Dufrêne}{Dufrene}, M. and P. Legendre. 1997. Species assemblages and indicator species: The need for a flexible asymetrical approach. Ecological Monographs 67:345-366.
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CTFC

Florian Jansen, Institute of Botany and Landscape Ecology, Ernst-Moritz-Arndt-University
}
\seealso{\code{\link{summary.multipatt}}, \code{\link{strassoc}}, \code{\link{signassoc}}, \code{\link{how}}}
\examples{
library(stats)

data(wetland) ## Loads species data

wetkm = kmeans(wetland, centers=3) ## Creates three clusters using kmeans

## Runs the combination analysis using IndVal.g as statistic
wetpt = multipatt(wetland, wetkm$cluster, control = how(nperm=999)) 

## Lists those species with significant association to one combination
summary(wetpt) 

## Lists those species with significant association to one combination, 
## including indval components.
summary(wetpt, indvalcomp=TRUE) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.