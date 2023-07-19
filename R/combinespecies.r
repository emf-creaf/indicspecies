#' Combines species from a community table
#' 
#' Creates species combinations to be used in indicator value analyses. 
#'
#' @param X A community data table with sites in rows and species in columns. This table can contain either presence-absence or abundance data.
#' @param min.order Minimum number of species conforming species combinations.
#' @param max.order Maximum number of species conforming species combinations.
#' @param min.occ Threshold for minimum occurrence. Combinations with lower values are not kept.
#' @param FUN Function to be calculated for all species forming the species combination. Accepted values are \code{min}, \code{max}, \code{sum} and \code{mean}.
#' @param verbose If TRUE, prints the results of each step.
#' @param add.names If TRUE, adds the names of the species combinations to the output matrix. Species combination names are lists of species concatenated using character '+'.
#' @param ... Additional arguments for function \code{FUN}.
#' 
#' @details
#' This function allows creating a data table where rows are sites and columns are combinations of species. Values for a given column of this matrix are derived from the abundance values of the species forming the corresponding combination. In particular, the abundance value for a given combination in a given site is equal to the value returned by function 'FUN' (normally the minimum) among the site values of all species forming the combination. The matrix 'XC' returned by this function can be used in functions \code{\link{strassoc}} and \code{\link{signassoc}}. Alternatively, \code{\link{indicators}} and related functions provide a more elaborated way to explore the indicator value of the simultaneous occurrence of sets of species (i.e. species combinations).
#' 
#'
#' @return
#' An list with:
#' \itemize{
#'   \item{XC - A matrix containing the abundance/occurrence of each species combination.}
#'   \item{C - A binary matrix describing the set of species forming each combination ('0' means that the species is not included, and '1' means that the species is included).}
#' }
#' 
#' @references 
#' De \enc{Cáceres}{Caceres}, M., Legendre, P., Wiser, S.K. and Brotons, L. 2012. Using species combinations in indicator analyses. Methods in Ecology and Evolution 3(6): 973-982.
#' 
#' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, EMF-CREAF
#' 
#' @seealso \code{\link{indicators}}, \code{\link{strassoc}}, \code{\link{signassoc}}
#' 
#' 
#' @export
#'
#' @examples
#' library(stats)
#' 
#' ## Loads species data
#' data(wetland) 
#' 
#' ## Create species combinations
#' Y=combinespecies(X=wetland, max.order=3, min.occ=5, verbose=TRUE)
#' 
#' ## Creates three site groups using kmeans
#' wetkm = kmeans(wetland, centers=3) 
#' 
#' ## Calculate indicator value of species combinations for each of the three site groups
#' strassoc(Y$XC, cluster=wetkm$cluster,func="IndVal.g") 
#' 
#' ## Calculate point biserial correlation value of species combinations 
#' ## for each of the three site groups
#' strassoc(Y$XC, cluster=wetkm$cluster,func="r.g") 
combinespecies<-function(X, min.order = 1, max.order = 3, min.occ = 1, FUN = min, verbose=FALSE, add.names=TRUE, ...) {
  ## Turn into a matrix (if not)
  X = as.matrix(X)
  
  if(sum(is.na(X))>0) stop("Cannot deal with NA values. Remove and run again.")
  if(min.order<1) stop("min.order must be at least 1.")
  if(min.order>max.order) stop("min.order must be equal or larger than max.order.")
  nsites = nrow(X)
  
  #Get species names
  spplist = colnames(X)
  nspecies = length(spplist)
  if(verbose) cat(paste("Number of species to combine: ",nspecies,"\n", sep=""))
  if(length(spplist)==1) stop("At least two species are necessary.")
  
  evalCombs<-function(spvec, dvec, verbose=FALSE) {
    comblist<-vector("list",0)
    sc.ab<-apply(X[,spvec, drop=FALSE],1,FUN=FUN,...)
    if(sum(ifelse(sc.ab>0,1,0))>=min.occ) {
      if(length(spvec)>=min.order) comblist<-c(comblist,list(spvec))
      if(length(spvec)<max.order && length(dvec)>0) { # explore other
        for(j in 1:length(dvec)) {
          if(verbose) cat(paste("Starting species ",dvec[j],"..."))
          comblist<-c(comblist,evalCombs(c(spvec,dvec[j]), dvec[-(1:j)], verbose=FALSE))
          if(verbose) cat(paste(" accepted combinations:",length(comblist),"\n"))
        }
      }
    }
    return(comblist)
  }
  
  veck = (1:length(spplist))#Vector of species indices to combine
  
  comblistDef<-vector("list",0)
  for(i in 1:length(veck)) {
    if(verbose) cat(paste("Starting species ",veck[i],"..."))
    comblistDef<-c(comblistDef,evalCombs(veck[i], veck[-(1:i)], verbose=FALSE))
    if(verbose) cat(paste(" accepted combinations:",length(comblistDef),"\n"))
  }
  nc = length(comblistDef)
  if(verbose) cat(paste("Final number of combinations: ",nc,"\n", sep=""))
  
  C = matrix(0,nrow=nc,ncol=nspecies)
  colnames(C)<-spplist
  XC = matrix(0, nrow=nrow(X), ncol=nc)
  rownames(XC) = rownames(X)
  for(i in 1:nc) {
    spvec = as.numeric(comblistDef[[i]])
    XC[,i] <-apply(X[,spvec, drop=FALSE],1,FUN=FUN,...)
    C[i,spvec] = 1
  }
  cn = 1:ncol(XC)
  if(add.names) {
    for(i in 1:nc) cn[i] = paste(spplist[C[i,]==1], collapse="+")
  }
  colnames(XC) = cn
  rownames(C) = cn
  return(list(XC = XC, C=C))
}