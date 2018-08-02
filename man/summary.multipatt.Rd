\encoding{UTF-8}
\name{summary.multipatt}
\alias{summary.multipatt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summarizing multi-level pattern analysis}
\description{
 This is a \code{summary} function for \code{multipatt} objects. 
}
\usage{
\method{summary}{multipatt}(object, alpha=0.05, minstat=NULL, At = NULL, Bt=NULL, indvalcomp=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object returned by function \code{\link{multipatt}}.}
  \item{alpha}{Significance level for selecting species in the summary.}
  \item{minstat}{Minimum value of the statistic for selecting species in the summary.}
  \item{At}{Minimum value of positive predictive value (A) for selecting indicators. This argument is effective only if \code{multipatt} function was executed using \code{func = "IndVal"} or \code{func = "IndVal.g"}.}
  \item{Bt}{Minimum value for sensitivity (B) for selecting indicators. This argument is effective only if \code{multipatt} function was executed using \code{func = "IndVal"} or \code{func = "IndVal.g"}.}
  \item{indvalcomp}{Flag used to report components "A" (or "A.g") and "B" of indicator values. This argument is effective only if \code{multipatt} function was executed using \code{func = "IndVal"} or \code{func = "IndVal.g"}.}
  \item{...}{Additional arguments of the \code{summary} function.}
}
\details{
   The summary function for \code{multipatt} objects presents the results of the analysis by grouping species following their associated site group combination. Only significantly associated species are listed for each combination. The user can limit the species shown by setting a significance level and the minimum value for the strength of the association. 
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Biodiversity and Landscape Ecology Laboratory, Forest Science Center of Catalonia
}
\seealso{\code{\link{strassoc}}, \code{\link{signassoc}}, \code{\link{multipatt}}}
\examples{
library(stats)

## Loads species data
data(wetland) 

## Creates three clusters using kmeans
wetkm = kmeans(wetland, centers=3) 

## Runs the combination analysis using IndVal.g as statistic
wetpt = multipatt(wetland, wetkm$cluster) 

## Lists those species with significant association to one combination
summary(wetpt) 

## Lists those species with significant association to one combination, 
## including indval components
summary(wetpt, indvalcomp=TRUE) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.