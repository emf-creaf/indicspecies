\encoding{UTF-8}
\name{signassoc}
\alias{signassoc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Statistical significance of species site-group associations}
\description{
 This function computes the permutation p-value of the association between a species vector and a vector of memberships to a site-group. Four different tests of hypotheses arise from considering either presence/absence or quantitative species data, and from using a non-equalized or group-equalized approach.
}
\usage{
signassoc(X, U = NULL, cluster = NULL, mode = 1, alternative = "greater", 
          control = how(), print.perm = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{species data (rows are sites, columns are species).}
  \item{U}{A matrix of group memberships. Each column corresponds to one site-group. Membership values can be crisp or fuzzy. If this matrix is not provided, vector gmv will be used.}
  \item{cluster}{a vector of numeric group labels for the sites.}
  \item{mode}{0 (site-based) or 1 (group-based).}
  \item{alternative}{Alternative hypothesis: "greater" (positive association), "less" (negative association), "two.sided" (either positive or negative).}
  \item{control}{a list of control values describing properties of the permutation design, as returned by a call to \code{\link{how}}.}
  \item{print.perm}{If TRUE, prints permutation numbers after each set of 100 permutations.}
}
\details{
Input data for this function is the species matrix X and either a matrix of group memberships (U) or a vector of numeric group labels (cluster). This R function works for both presence/absence and quantitative species data, depending on the values of the input matrix X.
If \code{mode = 0}, the null ecological hypothesis is that the frequency (or abundance) of the species of interest in sites belonging to the site group is not higher than the frequency (or abundance) in sites not belonging to it. If \code{mode = 1}, the null ecological hypothesis is that the relative frequency (average abundance) of the species of interest is not higher in the target site group than in other groups. See De Cáceres and Legendre for more details. Complex permutation designs are allowed through the function \code{\link{how}} from package "permute".
}
\value{
Returns a matrix of p-values, where species are in rows and groups are in columns. Two additional columns indicate the group with lowest p-value and the p-value for this group after Sidak's correction for multiple testing. 
}
\note{
Users should be aware that the significance test in \code{signassoc} is not exactly the same as the one in \code{indval} from \code{labdsv} package. The \code{signassoc} function is using the preference for the target group (either non-equalized or group-equalized) as test statistic. After every permutation the preference for the target group is recalculated. The function is therefore testing the null hypothesis stating that the preference of the species for a given site group is due to chance only (as in Bakker 2008). The test is repeated for every group, and this is the reason why there are as many p-values as groups. In contrast, the \code{indval} function from \code{labdsv} package uses the maximum preference value as test statistic, and the maximum preference value is recalculated after each permutation. The maximum preference may correspond to other groups than the one chosen for the unpermuted data. \code{indval} function from \code{labdsv} package is therefore testing the null hypothesis saying that the group with observed maximum preference is not such, because the maximum preference was in that group due to chance only. In order to get the consistent results compared to the \code{indval} function, users should use the function \code{\link{multipatt}}, along with the option \code{duleg=TRUE}. 
}
\references{
Bakker, J. 2008. Increasing the utility of Indicator Species Analysis. Journal of Applied Ecology 45: 1829-1835.

De \enc{Cáceres}{Caceres}, M. and Legendre, P. 2009. Associations between species and groups of sites: indices and statistical inference. Ecology 90(12): 3566-3574.
}

\author{Miquel De \enc{Cáceres}{Caceres} Ainsa,  CTFC}

\seealso{\code{\link{strassoc}}, \code{\link{multipatt}}, \code{\link{how}}}
\examples{
library(stats)

data(wetland) ## Load species data

wetkm = kmeans(wetland, centers=3)## Create three clusters using

## Look for species whose abundance is significantly higher in one of the three groups
signassoc(wetland, cluster=wetkm$cluster, mode=1, control = how(nperm=999))

## Look for species whose abundance is significantly higher in sites belonging 
## to one group as opposed to sites not belonging to it.
signassoc(wetland, cluster=wetkm$cluster, mode=0, control = how(nperm=999)) 
}
