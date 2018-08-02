\encoding{UTF-8}
\name{coverage}
\alias{coverage}
\alias{plotcoverage}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Coverage of a set of indicators}
\description{
Function \code{coverage} calculates the proportion of sites of the target site group where one or another indicator (a species or a species combination) is found. Parameters are used to select valid indicators before calculating coverage. Function \code{plotcoverage} plots the coverage against the positive predictive value threshold used to select indicators.
}
\usage{
coverage(x, y=NULL, selection=NULL, minstat=NULL, At=NULL, Bt=NULL, type="stat", 
alpha=0.05) 
plotcoverage(x, y=NULL, by=0.05, type="stat", max.order=NULL, group = NULL,  
alpha = 0.05, add=FALSE, xlab=expression(A[t]), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class '\code{\link{indicators}}' or a community data table with sites in rows and species in columns (in this second case, an object of class '\code{\link{multipatt}}' must be supplied for \code{y}).}
  \item{y}{An object of class '\code{\link{multipatt}}'.}
  \item{selection}{A logical vector restricting the set of indicators used to calculate the coverage.}
  \item{minstat}{Minimum value of the statistic for selecting indicators.}
  \item{At}{Minimum value of positive predictive value (A) for selecting indicators.}
  \item{Bt}{Minimum value for sensitivity (B) for selecting indicators.}
  \item{alpha}{Significance level for selecting indicators. This parameter is only used when the function is called using objects of class '\code{multipatt}'.}
  \item{type}{Specifies how to select indicators: (1) using the value of the statistic (\code{type = "stat"}); (2) the lower bound of its confidence interval (\code{type = "lowerCI"}); or (3) the upper bound of its confidence interval (\code{type = "upperCI"}). This parameter makes sense when the function is called using objects of class '\code{indicators}' and bootstrap confidence intervals are available for this object. Otherwise \code{type} has no effect and the value of the statistic is used for selection. In function \code{coverage}, the value of \code{type} applies to selection using \code{minstat}, \code{At} and \code{Bt}. In function \code{plotcoverage}, the value of \code{type} applies to selection using \code{At}.}
  \item{by}{Rate of increase in the predictive value threshold (At).}
  \item{max.order}{The maximum number of species conforming species combinations (allows examining the effects of increasing the order of combinations).This parameter is only used when the function is called using objects of class '\code{indicators}'.}
  \item{group}{Either an integer or a character string indicating the site group or site group combination for which plot is desired. This parameter is only used when the function is called using objects of class '\code{multipatt}'.}
  \item{add}{Flag indicating whether the function should draw on the existing plot.}
  \item{xlab}{Label for the x-axis.}
  \item{...}{Additional plotting parameters that are passed to the \code{plot} function.}
}
\details{
The \code{coverage} of a set of indicators was defined in De \enc{Cáceres}{Caceres} et al. (2012) as the proportion of sites in a given site group where one or several indicators are found. This value allows assessing how often the site group will be able to be determined. If all indicators of a site group are rare, then the indication system will not be useful, regardless of how much restricted to the site group the indicators are. The coverage value is a generalization of quantity B of IndVal, that applies to a group of indicators instead of a single one. Function \code{plotcoverage} plots the coverage against the positive predictive value threshold (At) used to select indicators, as in De \enc{Cáceres}{Caceres} et al. (2012). Functions \code{coverage} and \code{plotcoverage} can be executed using either an object of class '\code{indicators}', or an object of class '\code{multipatt}'. However, the parameters that apply to each case are slightly different. When using \code{coverage} and \code{plotcoverage} on objects of class '\code{multipatt}' one is expected to calculate the coverage for those indicators that are significant (see \code{alpha} parameter), although other constraints to select valid indicators can be used. When using \code{coverage} and \code{plotcoverage} on objects of class '\code{indicators}' one is expected to calculate the coverage for indicators that have values of A larger than a specified threshold (see \code{At} parameter). In this latter case, it may be advisable to use \code{stat="lowerCI"}, so that indicators with broad confidence intervals are not included in the selection.
}
\value{
  When used with an object of class '\code{indicators}', function \code{coverage} returns the proportion of sites of the target site group where one or another indicator (species combination) is found. When used with an object of class '\code{indicators}', function \code{coverage} returns a vector containing the coverage value for each site group or site group combination.
}
\references{
De \enc{Cáceres}{Caceres}, M., Legendre, P., Wiser, S.K. and Brotons, L. 2012. Using species combinations in indicator analyses. Methods in Ecology and Evolution 3(6): 973-982.
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Biodiversity and Landscape Ecology Laboratory, Centre Tecnologic Forestal de Catalunya
}
\seealso{\code{\link{indicators}}, \code{\link{multipatt}}, \code{\link{pruneindicators}}}
\examples{
library(stats)

data(wetland) ## Loads species data

## Creates three clusters using kmeans
wetkm = kmeans(wetland, centers=3) 

## Determine sensitivity of individual species
B=strassoc(wetland, cluster=wetkm$cluster,func="B") 

## Select species with more than 20% of sensitivity for the first group
sel=which(B[,1]>0.2) 

## Run indicator analysis with species combinations for the first group
sc= indicators(X=wetland[,sel], cluster=wetkm$cluster, group=1, verbose=TRUE, At=0.5, Bt=0.2)

## Determine the coverage of the selected set of indicators
coverage(sc)

## Plot the coverage against the threshold At
plotcoverage(sc)
plotcoverage(sc, max.order=2, add=TRUE, lty=2)

## Runs the combination analysis using IndVal.g as statistic
wetpt = multipatt(wetland, wetkm$cluster, control = how(nperm=999)) 
## Determines the coverage for each site group combination
coverage(wetland, wetpt)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.