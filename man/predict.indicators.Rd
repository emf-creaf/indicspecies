\encoding{UTF-8}
\name{predict.indicators}
\alias{predict.indicators}
\alias{predict.indicators.cv}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Predicts site group from indicators}
\description{
 Function \code{predict.indicators} takes an object \code{\link{indicators}} and determines the probability of the indicated site group given a new community data set. Function \code{predict.indicators.cv} allows estimating the cross-validated (by leave-one-out) probability of each site belonging to the target site group.
}
\usage{
\method{predict}{indicators} (object, newdata,...)
predict.indicators.cv(X, cluster, group, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class 'indicators'.}
  \item{newdata}{A community data table (with sites in rows and species in columns) for which predictions are needed. This table can contain either presence-absence or abundance data, but only presence-absence information is used for the prediction.}
  \item{...}{In function \code{predict}, additional arguments not used (included for compatibility with \code{\link{predict}}). In function \code{predict.indicators.cv}, additional arguments passed onto function \code{\link{indicators}}.}
  \item{X}{A community data table with sites in rows and species in columns. This table can contain either presence-absence or abundance data.}
  \item{cluster}{A vector containing the classification of sites into site groups.}
  \item{group}{The label corresponding to the target site group. }
}
\details{
   Function \code{\link{indicators}} explores the indicator value of the simultaneous occurrence of sets of species (i.e. species combinations). The method is described in De \enc{Cáceres}{Caceres} et al. (2012) and is a generalization of the Indicator Value method of \enc{Dufrêne}{Dufrene} & Legendre (1997). The current function \code{predict.indicators} is used to predict the indicated site group from the composition of a new set of observations. For communities where one or more of the indicator species combinations are found, the function returns the probability associated to the indicator that has the highest positive predictive value (if confidence intervals are available, the maximum value is calculated across the lower bounds of the confidence interval). For communities where none of the indicator species combinations is found, the function returns zeroes. 
   
   Function \code{predict.indicators.cv} is used to evaluate the predictive power of a set of indicators in a cross-validated fashion. For each site in the data set, the function makes a call to function \code{\link{indicators}} excluding the information of the site, and then calls \code{predict.indicators} to evaluate the probability of the site group. 
}
\value{
  If confidence intervals are available in \code{x}, function \code{predict.indicators} returns a matrix where communities are in rows and there are three columns, correspoinding to the probability of the indicated site group along with the confidence interval. If confidence intervals are not available in \code{x}, then \code{predict.indicators} returns a single vector with the probability of the indicated site group for each community. Function \code{predict.indicators.cv} returns a vector with the (cross-validated) probability of the indicated site group for each community. 
}
\references{
De \enc{Cáceres}{Caceres}, M., Legendre, P., Wiser, S.K. and Brotons, L. 2012. Using species combinations in indicator analyses. Methods in Ecology and Evolution 3(6): 973-982.

\enc{Dufrêne}{Dufrene}, M. and P. Legendre. 1997. Species assemblages and indicator species: The need for a flexible asymetrical approach. Ecological Monographs 67:345-366.

}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Biodiversity and Landscape Ecology Laboratory, Forest Science Center of Catalonia
}
\seealso{\code{\link{indicators}}, \code{\link{pruneindicators}} \code{\link{coverage}}, \code{\link{multipatt}}, \code{\link{strassoc}}, \code{\link{signassoc}}}
\examples{
library(stats)

data(wetland) ## Loads species data

## Creates three clusters using kmeans
wetkm = kmeans(wetland, centers=3) 

## Determine sensitivity of individual species
B=strassoc(wetland, cluster=wetkm$cluster,func="B") 

## Select species with more than 20\% of sensitivity for the first group
sel=which(B[,1]>0.2) 

## Run indicator analysis with species combinations for the first group
sc= indicators(X=wetland[,sel], cluster=wetkm$cluster, group=1, verbose=TRUE, At=0.5, Bt=0.2)

## Use the indicators to make predictions of the probability of group #1
## Normally an independent data set should be used, because 'wetland' was used to derive
## indicators.
p<-predict(sc, wetland)

## Calculate cross-validated probabilities (calls functions indicators once for each site 
## after excluding it, and then calls predict.indicators for that site)
pcv = predict.indicators.cv(X=wetland[,sel], cluster=wetkm$cluster, group=1, At=0.5, Bt=0.2)

## Show original membership to group 1 along with (resubstitution) predicted probabilities  
## and cross-validated probabilities
data.frame(Group1 = as.numeric(wetkm$cluster==1), Prob = p, Prob_CV = pcv)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.