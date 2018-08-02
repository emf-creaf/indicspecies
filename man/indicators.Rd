\encoding{UTF-8}
\name{indicators}
\alias{indicators}
\alias{print.indicators}
\alias{summary.indicators}
\alias{plot.indicators}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Indicator analysis for species combinations}
\description{
 Determines the indicator value of species combinations. 
}
\usage{
indicators(X, cluster, group, func="IndVal", max.order=5, max.indicators=NULL, 
At=0, Bt=0, sqrtIVt=0, nboot=0, alpha=0.05, XC=TRUE, enableFixed = FALSE, verbose = FALSE)
\method{print}{indicators} (x, At=0, Bt=0, sqrtIVt=0, selection=NULL, confint=FALSE,...)
\method{plot}{indicators} (x, type="sqrtIV", maxline=TRUE,...)
\method{summary}{indicators} (object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{A community data table with sites in rows and species in columns. This table can contain either presence-absence or abundance data.}
  \item{x}{An object of class 'indicators'.}
  \item{object}{An object of class 'indicators'.}
  \item{cluster}{A vector containing the classification of sites into site groups.}
  \item{group}{The label corresponding to the target site group. }
  \item{max.order}{Maximum number of species conforming species combinations.}
  \item{max.indicators}{Maximum number of valid indicators to be kept. If \code{NULL}, then all valid indicators are kept.}
  \item{func}{The indicator value variant to be used, either "IndVal" (non-equalized) or "IndVal.g" (group-equalized).}
  \item{At}{Threshold for positive predictive value used to select valid indicators. Combinations with lower values are not kept.}
  \item{Bt}{Threshold for sensitivity used to select valid indicators. Combinations with lower values are not kept.}
  \item{sqrtIVt}{Threshold for (square root of) indicator value. Combinations with lower values are not kept.}
  \item{nboot}{Number of bootstrap samples for confidence intervals.}
  \item{alpha}{Error in confidence intervals.}
  \item{XC}{If TRUE, outputs the abundance/occurrence matrix of species combinations.}
  \item{enableFixed}{If TRUE, uses species that occur in all sites as fixed elements and creates combinations with the remaining ones.}
  \item{verbose}{If TRUE, prints the results of each step.}
  \item{selection}{A logical vector used to restrict, a priori, the species combinations to be printed.}
  \item{confint}{Flag to indicate that confidence interval bounds are desired when printing.}
  \item{type}{Statistic to plot. Accepted values are "IV" (indicator value), "sqrtIV" (square root of indicator value), "A", "LA", "UA", (positive predictive value and confidence limits), "B",  "LB", "UB" (sensitivity and confidence limits).}
  \item{maxline}{Flag to indicate whether a line has to be drawn joining the maximum values for each order of combinations.}
  \item{...}{Additional arguments for functions \code{print},\code{summary} or \code{plot}.}
}
\details{
   Function \code{indicators} creates explores the indicator value of the simultaneous occurrence of sets of species (i.e. species combinations). The method is described in De \enc{Cáceres}{Caceres} et al. (2012) and is a generalization of the Indicator Value method of \enc{Dufrêne}{Dufrene} & Legendre (1997). The maximum number of species conforming the species combination can be controlled using \code{max.order}. For each combination of species it determines its positive predictive value (A), sensitivity (B) and the square root of indicator value (sqrtIV). If \code{nboot>0} then bootstrap confidence intervals are determined with the specified \code{alpha} level, as explained in De \enc{Cáceres}{Caceres} & Legendre (2009). The combinations to be kept can be restricted to those whose positive predictive value, sensitivity and/or indicator value are equal or greater than input thresholds. Function \code{print} allows printing the results in a nice table, whereas \code{summary} provides information about candidate species, combinations and coverage of the set of indicators. Function \code{plot} draws the statistics against the order (i.e. the number of species) of the combination.
}
\value{
  An object of class \code{indicators} with:
  \item{candidates}{The vector of initial candidate species.}
  \item{finalsplist}{The vector of species finally selected for combinations.}
  \item{C}{A matrix describing all the combinations studied.}
  \item{XC}{A matrix containing the abundance/occurrence of each species combination.}
  \item{A}{Positive predictive power of species combinations. If \code{nboot>0} then this includes the lower and upper bounds of the confidence interval.}
  \item{B}{Sensitivity of species combinations. If \code{nboot>0} then this includes the lower and upper bounds of the confidence interval.}
  \item{sqrtIV}{Square root of indicator value of species combinations. If \code{nboot>0} then this includes the lower and upper bounds of the confidence interval.}
  \item{group.vec}{A logical vector indicating the membership to the target group.}
}
\references{
De \enc{Cáceres}{Caceres}, M., Legendre, P., Wiser, S.K. and Brotons, L. 2012. Using species combinations in indicator analyses. Methods in Ecology and Evolution 3(6): 973-982.

De \enc{Cáceres}{Caceres}, M. and Legendre, P. 2009. Associations between species and groups of sites: indices and statistical inference. Ecology 90(12): 3566-3574.

\enc{Dufrêne}{Dufrene}, M. and P. Legendre. 1997. Species assemblages and indicator species: The need for a flexible asymetrical approach. Ecological Monographs 67:345-366.

}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Biodiversity and Landscape Ecology Laboratory, Forest Science Center of Catalonia
}
\seealso{\code{\link{predict.indicators}},\code{\link{pruneindicators}}, \code{\link{coverage}}, \code{\link{multipatt}}, \code{\link{strassoc}}, \code{\link{signassoc}}}
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
sc= indicators(X=wetland[,sel], cluster=wetkm$cluster, group=1, verbose=TRUE, 
               At=0.5, Bt=0.2)

#Prints the results
print(sc)

## Plots positive predictive power and sensitivity against the order of 
## combinations
plot(sc, type="A")
plot(sc, type="B")

## Run indicator analysis with species combinations for the first group, 
## but forcing 'Orysp' to be in all combinations
sc2= indicators(X=wetland[,sel], cluster=wetkm$cluster, group=1, verbose=TRUE, 
                At=0.5, Bt=0.2, enableFixed=TRUE)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.