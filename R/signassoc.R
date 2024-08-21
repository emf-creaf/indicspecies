#' Statistical significance of species site-group associations
#' 
#' This function computes the permutation p-value of the association between a species vector and a vector of memberships to a site-group. Four different tests of hypotheses arise from considering either presence/absence or quantitative species data, and from using a non-equalized or group-equalized approach.
#'
#' @param X Community data table (rows are sites, columns are species).
#' @param U A matrix of group memberships. Each column corresponds to one site-group. Membership values can be crisp or fuzzy. If this matrix is not provided, vector gmv will be used.
#' @param cluster A vector of numeric group labels for sites.
#' @param mode Association mode, 0 (site-based) or 1 (group-based).
#' @param alternative Alternative statistical hypothesis: "greater" (positive association), "less" (negative association), "two.sided" (either positive or negative).
#' @param control A list of control values describing properties of the permutation design, as returned by a call to \code{\link[permute]{how}}.
#' @param permutations A custom matrix of permutations, to be used if \code{control = NULL}, with permutations in rows and site indices in columns.
#' @param print.perm If TRUE, prints permutation numbers after each set of 100 permutations.
#'
#' @details
#' Input data for this function is the species matrix X and either a matrix of group memberships (U) or a vector of numeric group labels (cluster). This R function works for both presence/absence and quantitative species data, depending on the values of the input matrix X.
#' If \code{mode = 0}, the null ecological hypothesis is that the frequency (or abundance) of the species of interest in sites belonging to the site group is not higher than the frequency (or abundance) in sites not belonging to it. If \code{mode = 1}, the null ecological hypothesis is that the relative frequency (average abundance) of the species of interest is not higher in the target site group than in other groups. See De Cáceres and Legendre for more details. 
#' 
#' Complex permutation designs are allowed through the function \code{\link[permute]{how}} from package "permute". If those are not enough, the user can set \code{control = NULL} and specify a custom matrix of permutations to test with parameter \code{permutations}.
#' 
#' @return
#' Returns a matrix of p-values, where species are in rows and groups are in columns. Two additional columns indicate the group with lowest p-value and the p-value for this group after Sidak's correction for multiple testing. 
#' 
#' @note
#' Users should be aware that the significance test in \code{signassoc} is not exactly the same as the one in \code{indval} from \code{labdsv} package. The \code{signassoc} function is using the preference for the target group (either non-equalized or group-equalized) as test statistic. After every permutation the preference for the target group is recalculated. The function is therefore testing the null hypothesis stating that the preference of the species for a given site group is due to chance only (as in Bakker 2008). The test is repeated for every group, and this is the reason why there are as many p-values as groups. In contrast, the \code{indval} function from \code{labdsv} package uses the maximum preference value as test statistic, and the maximum preference value is recalculated after each permutation. The maximum preference may correspond to other groups than the one chosen for the unpermuted data. \code{indval} function from \code{labdsv} package is therefore testing the null hypothesis saying that the group with observed maximum preference is not such, because the maximum preference was in that group due to chance only. In order to get the consistent results compared to the \code{indval} function, users should use the function \code{\link{multipatt}}, along with the option \code{duleg=TRUE}.
#' 
#' @references 
#' Bakker, J. 2008. Increasing the utility of Indicator Species Analysis. Journal of Applied Ecology 45: 1829-1835.
#' 
#' De \enc{Cáceres}{Caceres}, M. and Legendre, P. 2009. Associations between species and groups of sites: indices and statistical inference. Ecology 90(12): 3566-3574.
#' 
#' @author Miquel De \enc{Cáceres}{Caceres} Ainsa,  EMF-CREAF
#' 
#' @seealso \code{\link{strassoc}}, \code{\link{multipatt}}
#' 
#' @export
#'
#' @examples
#' library(stats)
#' 
#' data(wetland) ## Load species data
#' 
#' wetkm <- kmeans(wetland, centers=3)## Create three clusters using
#' 
#' ## Look for species whose abundance is significantly higher in one of the three groups
#' signassoc(wetland, cluster=wetkm$cluster, mode=1, control = how(nperm=999))
#' 
#' ## Look for species whose abundance is significantly higher in sites belonging 
#' ## to one group as opposed to sites not belonging to it.
#' signassoc(wetland, cluster=wetkm$cluster, mode=0, control = how(nperm=999)) 
#' 
signassoc <- function(X, U=NULL, cluster=NULL, mode = 1, alternative="greater", 
                      control = how(), permutations = NULL, print.perm=FALSE) {
	
  vector.to.partition <- function(v, clnames) {
    m <- t(sapply(v,function(x) as.numeric(x==clnames)))
    dimnames(m) = list(1:length(v),clnames)
    return(m)                                                                                                                              
  }                                                                                                                                          
  

  #Turn into a matrix (if not)
  X = as.matrix(X)

  nsps = ncol(X)
  spnames = colnames(X)
  nsites = nrow(X)

  if(!is.null(control)){
    nperm = control$nperm
  } else if (!is.null(permutations)){
    nperm = nrow(permutations)
  } else {
    stop("You must control permutations with how() or supply a matrix of permutations")
  }
  
   mode= match.arg(as.character(mode), c("0","1"))
   alternative= match.arg(as.character(alternative), c("greater","less","two.sided"))
  if(sum(is.na(X))>0) stop("Cannot deal with NA values. Remove and run again.")
  
  if(is.null(U)) U = vector.to.partition(cluster, levels(factor(cluster))) 
  if(sum(is.na(U))>0) stop("Cannot deal with NA values. Remove and run again.")
  
  ngroups = ncol(U)	
  cdm = matrix(1,nrow=nsps,ncol=ngroups)
  ddm = matrix(1,nrow=nsps,ncol=ngroups)
  
  U = as.matrix(U)
  if(mode==0) {
	  dm = t(X)%*%U
  } else {
	  aisp = t(X)%*%U
  	  ni = diag(t(U)%*%U)
  	  aispni=sweep(aisp,2,ni,"/")
  	  aispni[is.na(aispni)]=0 # check for division by zero
  	  s = apply(aispni,1,"sum")
  	  dm = sweep(aispni,1,s,"/")
  	  dm[is.na(dm)]=0 # check for division by zero
  }
  	
  a <- system.time({
  for(p in 1:nperm) {
      if(p%%100==0 & print.perm) cat("perm", p,"\n")
    if(!is.null(control)){
      pInd = shuffle(nsites, control=control)
    } else {
      pInd = permutations[p,]
    }
		pX = as.matrix(X[pInd,])
  		if(mode==0) {
	  		dmp = t(pX)%*%U
  		} else {
	 		aisp = t(pX)%*%U
		  	ni = diag(t(U)%*%U)
  	  		aispni=sweep(aisp,2,ni,"/")
	    	aispni[is.na(aispni)]=0 # check for division by zero
  	  		s = apply(aispni,1,"sum")
  	  		dmp = sweep(aispni,1,s,"/")
	   	   dmp[is.na(dmp)]=0 # check for division by zero
  		}			
		if(alternative=="less") {
		   cdm = cdm + as.numeric(dmp<=dm)
		} else if(alternative=="greater") {
		   cdm = cdm + as.numeric(dmp>=dm)
		} else if(alternative=="two.sided") {
		   cdm = cdm + as.numeric(dmp>=dm)
		   ddm = ddm + as.numeric(dmp<=dm)
		}
	}
	if(alternative!="two.sided") {
		cdm=cdm/(nperm+1)
	} else {
		cdm=pmin(matrix(1,nrow=nsps,ncol=ngroups),(2*pmin(cdm,ddm))/(nperm+1))
	}
	rownames(cdm)=spnames
  colnames(cdm)=colnames(U)
	})
  a[3] <- sprintf("%2f",a[3])
  #cat("Time to compute p-values =",a[3]," sec",'\n')


   psidak = vector(mode="numeric", length=nsps)
   best = vector(mode="numeric", length=nsps)
   for(i in 1:nsps) {
   	  best[i] = which(cdm[i,]==min(cdm[i,]))[1]
	  psidak[i] = (1-(1-min(cdm[i,]))^ngroups)
   	}
	cdm = cbind(cdm,best,psidak)
	colnames(cdm) = c(levels(as.factor(cluster)),"best","psidak")
	return(cdm)      
 }

