\encoding{UTF-8}
\name{indicspecies-package}
\alias{indicspecies-package}
\alias{indicspecies}
\docType{package}
\title{
 Studying the statistical relationship between species and groups of sites 
}
\description{
This package provides a set of functions to assess the strength and statistical significance of the relationship between species occurrence/abundance and groups of sites. It is also possible to check the statistical significance of such associations.
New functions have been added (previously distributed in 'resniche') which allow evaluate several niche metrics while accounting for distances between resources on the basis of the similar perception from species (De \enc{Cáceres}{Caceres} et al. 2011).
}
\details{
\tabular{ll}{
Package: \tab indicspecies\cr
Type: \tab Package\cr
Version: \tab 1.7.7\cr
Date: \tab 2019-05-20\cr
Depends: \tab permute \cr
Suggests: \tab sp, rgeos \cr
License: \tab GPL (>=2)\cr
LazyLoad: \tab yes\cr
}
Use \code{\link{strassoc}} for assessing strength of associations and \code{\link{signassoc}} to assess the statistical significance of the associations (De \enc{Cáceres}{Caceres} & Legendre 2009). Much more complete functions are \code{\link{multipatt}} (De \enc{Cáceres}{Caceres} et al. 2010) and \code{\link{indicators}} (De \enc{Cáceres}{Caceres} et al. 2012). Read the package vignette using  \code{vignette("indicspeciesTutorial", package="indicspecies")} to learn how to use this package.

}
\author{
Miquel De \enc{Cáceres}{Caceres}, Florian Jansen

Maintainer: Miquel De \enc{Cáceres}{Caceres} <miquelcaceres@gmail.com>
}
\references{
De \enc{Cáceres}{Caceres}, M. and Legendre, P. 2009. Associations between species and groups of sites: indices and statistical inference. Ecology 90(12): 3566-3574.

De \enc{Cáceres}{Caceres}, M., Legendre, P., Moretti, M. 2010. Improving indicator species analysis by combining groups of sites. Oikos 119(10): 1674-1684.

De \enc{Cáceres}{Caceres}, M., Legendre, P., Wiser, S.K. and Brotons, L. 2012. Using species combinations in indicator analyses. Methods in Ecology and Evolution 3(6): 973-982.

De \enc{Cáceres}{Caceres}, M., Sol, D., Lapiedra, O. and P. Legendre. (2011) A framework for estimating niche metrics using the resemblance between qualitative resources. Oikos 120: 1341-1350.
}
\keyword{ package }