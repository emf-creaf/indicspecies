\name{nichevar}
\Rdversion{1.1}
\alias{nichevar}
\alias{nichecentroid}
\alias{nichepref}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Resource niche metrics for a single niche
}
\description{
Function \code{nichepref} computes the species resource preference from a the species resource use (and resource availability when given). Function \code{nichecentroid} computes the centroid on the resource space for a set of species. Function \code{nichevar} computes the multivariate resource variance for a set of species (i.e. niche breadth). In all functions resources are given in distance matrix \code{D}, the species resource use is given in \code{P} and the availability of resources, if present, are given by vector \code{q}. 
}
\usage{
nichevar(P, D = NULL, q = NULL, mode="multiple", Np = NULL, 
Nq = NULL, nboot = 1000, alpha=0.05)
nichecentroid(P, D = NULL, q = NULL, mode="multiple", Np = NULL, 
Nq = NULL, nboot = 1000, alpha=0.05)
nichepref(P, D = NULL, q = NULL, mode="multiple", Np = NULL, 
Nq = NULL, nboot = 1000, alpha=0.05)
}
\arguments{
  \item{P}{Data frame containing the relative or absolute usage that a set of species (in rows) make of a set of resources (in columns).}
  \item{D}{Object of type \code{\link{dist}} containing distance values between resources. If no distance matrix is provided (i.e. if \code{D==NULL}), the distances between resources is assumed to be maximum.}
  \item{q}{Vector with the availability of each resource.}
  \item{mode}{Either \code{mode = "single"} (rows of matrix P are individual observations to be pooled for a single niche) or \code{mode = "multiple"} (rows in P represent different niches).}
  \item{Np}{Vector with the number of observations per species from which the values in \code{P} come (in \code{mode = "multiple"}).}
  \item{Nq}{The number of observations per species from which the values in \code{q} come (in \code{mode = "multiple"}).}
  \item{nboot}{Number of boostrap samples used to compute bias-corrected percentile confidence intervals.}
  \item{alpha}{Used to set the confidence level (i.e. \code{alpha = 0.05} means 95 percent confidence interval).}
}

\details{
The method is described in De Caceres et al. (2010). If the distance matrix is not specified (i.e. if \code{D=NULL}) the function assumes that all resources are at a maximum distance (\code{d=1}). If the resource availability vector \code{q} is given then the values in \code{P} are taken as assessments of resource use and the species preference is calculated taking into account resource availability. Otherwise resource use is equated to resource preference. Moreover, most functions can compute bootstrap confidence intervals following the bias-corrected percentile method (Manly 2007). If \code{mode = "multiple"} and \code{Np != NULL}, bootstrap samples for a given species are generated assuming a multinomial distribution with the proportions calculated from the corresponding row values in \code{P}, and the number of observations comes from the corresponding element in \code{Np}. If \code{mode = "single"} then the bootstrapped units are the rows of matrix \code{P}. In both cases, if \code{Nq} is indicated the availability of resources is also bootstrapped. The bias-corrected percentile method was described for overlap niche measures in Mueller and Altenberg (1985) and is extended here for all niche metrics.
}

\value{
Function \code{nichepref} returns a matrix of species relative preference. Function \code{nichevar} returns a vector with the variance of the resources used for each species in \code{P}. Function \code{nichecentroid} returns a matrix niche centroid in the resource space for each species in \code{df}. If bootstrap confidence intervals are asked then the three functions also compute two extra data containing respectively the lower and upper bounds of the confidence intervals obtained following the bias-corrected percentile method. Function \code{nichearea} returns the area of the convex hull occupied by the resources used for each species in \code{P}. 
}

\references{
Mueller, L.D. and L. Altenberg. 1985. Statistical Inference on Measures of Niche Overlap. Ecology 66:1204-1210.

Manly, B.F.J. 2007. Randomization, bootstrap and Monte Carlo methods in biology. Chapman and Hall texts in statistical science series. 2nd edition.

De Caceres, M., Sol, D., Lapiedra, O. and P. Legendre. (2011) A framework for estimating niche metrics using the resemblance between qualitative resources. Oikos 120: 1341-1350.
}

\author{
Miquel De Caceres Ainsa, EMF-CREAF
}

\seealso{
See \code{\link{nicheoverlap}} for descriptors comparing two niches.
}
\examples{
# Loads example data
data(birds)

# The niche metrics using distances among resources and assuming equal availability of resources
nichepref(birdsbreed, D = resourceD) 
nichevar(birdsbreed, D = resourceD) 
nichecentroid(birdsbreed, D = resourceD) 

# The niche metrics using distances among resources and computes 
# 95 percent confidence intervals
nichepref(birdsbreed, D = resourceD, mode="multiple", 
Np = rowSums(birdsbreed), Nq = 100) 
nichevar(birdsbreed, D = resourceD, mode="multiple", 
Np = rowSums(birdsbreed), Nq = 100) 
nichecentroid(birdsbreed, D = resourceD, mode="multiple", 
Np = rowSums(birdsbreed), Nq = 100) 

# Same computations with different resource availability
nichepref(birdsbreed, D = resourceD, 
q = c(0.18, 0.24, 0.22, 0.21, 0.15), mode="multiple")
nichevar(birdsbreed, D = resourceD, 
q = c(0.18, 0.24, 0.22, 0.21, 0.15), mode="multiple")
nichecentroid(birdsbreed, D = resourceD, 
q = c(0.18, 0.24, 0.22, 0.21, 0.15), mode="multiple")

# The niche metrics using distances among resources and 
# computes 95 percent confidence intervals
nichepref(birdsbreed, D = resourceD, 
q = c(0.18, 0.24, 0.22, 0.21, 0.15), mode="multiple", Np = rowSums(birdsbreed), Nq = 100)
nichevar(birdsbreed, D = resourceD, 
q = c(0.18, 0.24, 0.22, 0.21, 0.15), mode="multiple", Np = rowSums(birdsbreed), Nq = 100)
nichecentroid(birdsbreed, D = resourceD, 
q = c(0.18, 0.24, 0.22, 0.21, 0.15), mode="multiple",  Np = rowSums(birdsbreed), Nq = 100)

}

