\encoding{UTF-8}
\name{strassoc}
\alias{strassoc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Strength of species site-group associations}
\description{
  This function computes the strength of the association between a species vector and a vector of memberships to the site-group. Several statistics are possible, following De \enc{Cáceres}{Caceres} and Legendre (2009). The choice of the most appropriate statistic depends on the application. If required, bootstrap confidence interval bounds are also computed.
}
\usage{
strassoc(X, cluster, func = "r", group = NULL, nboot = 0, alpha = 0.05, c = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{species data (rows are sites, columns are species).}
  \item{cluster}{a vector of numeric group labels for the sites.}
  \item{func}{the association index to be used. Accepted function names: "r", "r.g", "IndVal", "IndVal.g", "A", "A.g", "B", "cos", "cos.g", "r.ind", "r.ind.g", "s.ind", "s.ind.g" (lowercase values are also accepted).}
  \item{group}{Group for which association values are to be computed. If group=NULL, association values will be computed for all groups.}
  \item{nboot}{Number of bootstrap samples.}
  \item{alpha}{Error in confidence intervals }
  \item{c}{Total abundance per site (used in functions "r.ind", "r.ind.g", "s.ind", "s.ind.g" only).}
}
\details{
  This R function is applicable to both presence-absence and quantitative species data, depending on the values in the input matrix X.  Indices "r","r.g","cos" and "cos.g" correspond to the general correlation generalization ("r" is the point-biserial correlation coefficient), whereas "r.ind", "r.ind.g", "s.ind" and "s.ind.g" correspond to the individual-based generalization. Both approaches give the same results for binary (presence-absence) species data. Indicator value indices "IndVal" and "IndVal.g" are partial generalizations of the indices of the presence-absence indices. The "IndVal" index of \enc{Dufrêne}{Dufrene} & Legendre (1997) is called "IndVal.g" here and \code{strassoc} returns the square root of the original index. Indices "A","A.g" and "B" are the asymmetric quantities into which indicator values "IndVal" and "IndVal.g" can be decomposed. See De Cáceres and Legendre for details on the relationships between these indices and their usage. Bootstrap confidence intervals are computed using the simple percentile method (Manly 1997).
}
\value{
  Returns a matrix of association values, where species are in rows and groups are in columns. If 'nboot' is non zero, then a list is returned, where 'stat' contains the matrix of association values, and the lower and upper confidence limits are given in supplementary matrices named 'lowerCI' and 'upperCI' respectively.
}
\note{
This function gives the same association values as function \code{indval} in package "labdsv" when used setting \code{func="IndVal.g"}, excepting the fact that the square root IndVal values is returned instead of the original IndVal.
}
\references{
De \enc{Cáceres}{Caceres}, M. and Legendre, P. 2009. Associations between species and groups of sites: indices and statistical inference. Ecology 90(12): 3566-3574.

\enc{Dufrêne}{Dufrene}, M. and P. Legendre. 1997. Species assemblages and indicator species: The need for a flexible asymetrical approach. Ecological Monographs 67:345-366.

Manly, B. F. J. 1997. Randomization, bootstrap and Monte Carlo methods in biology. Chapman and Hall Texts in Statistical Science Series.
}
\author{Miquel De \enc{Cáceres}{Caceres} Ainsa,  Biodiversity and Landscape Ecology Laboratory, Centre Tecnologic Forestal de Catalunya}
\seealso{\code{\link{signassoc}}, \code{\link{multipatt}}}
\examples{
## Load species data
data(wetland) 

## Create three clusters using
wetkm = kmeans(wetland, centers=3)

## Compute Dufrene and Legendre's IndVal
strassoc(wetland, wetkm$cluster, func="IndVal.g") 

## Compute point-biserial correlation, with bootstrap 95 percent confidence intervals
strassoc(wetland, wetkm$cluster, func="r", nboot =100) 
}
