#' Determines the best subset of indicators
#' 
#' This function allows reducing drastically the number of species combinations to be retained for a given target site group.
#'
#' @param x An object of class '\code{\link{indicators}}'.
#' @param At Threshold for positive predictive value. Combinations with lower values are not kept.
#' @param Bt Threshold for sensitivity. Combinations with lower values are not kept.
#' @param sqrtIVt Threshold for (square root of) indicator value. Combinations with lower values are not kept.
#' @param alpha Threshold for statistical significance of indicator value. Combinations with higher p-values are not kept.
#' @param max.indicators Maximum number of species combinations to be kept. If \code{NULL}, the function returns all the non-nested valid indicators without further selection.
#' @param verbose If TRUE, prints the results of each step.
#' 
#' @details
#' First, the function selects those indicators (species or species combinations) with valid positive predictive value, sensitivity and indicator value, according to the input thresholds. If the object '\code{speciescomb}' contains confidence intervals, then the lower bounds are used to select the valid indicators. Second, the function discards those valid indicators whose occurrence pattern is nested within other valid indicators. Third, the function evaluates the \code{\link{coverage}} of the remaining set of indicators and explores subsets of increasing number of indicators, until the same coverage is attained and the set of indicators is returned. If the maximum allowed members is attained (\code{max.indicators}) then the set of indicators with maximum coverage is returned.
#' 
#'
#' @return
#' An object of class '\code{\link{indicators}}' with only the species combinations selected.
#' 
#' @references 
#' De \enc{Cáceres}{Caceres}, M., Legendre, P., Wiser, S.K. and Brotons, L. 2012. Using species combinations in indicator analyses. Methods in Ecology and Evolution 3(6): 973-982.
#' 
#' De \enc{Cáceres}{Caceres}, M. and Legendre, P. 2009. Associations between species and groups of sites: indices and statistical inference. Ecology 90(12): 3566-3574.
#' 
#' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, EMF-CREAF
#' 
#' @seealso \code{\link{indicators}}, \code{\link{coverage}}
#' 
#' 
#' @export
#'
#' @examples
#' library(stats)
#' 
#' data(wetland) ## Loads species data
#' 
#' ## Creates three clusters using kmeans
#' wetkm <- kmeans(wetland, centers=3) 
#' 
#' ## Run indicator analysis with species combinations for the first group
#' sc <- indicators(X=wetland, cluster=wetkm$cluster, group=1, verbose=TRUE, At=0.5, Bt=0.2)
#' 
#' ## Finds the 'best' subset of indicators
#' sc2 <- pruneindicators(sc, At=0.5, Bt=0.2, verbose=TRUE)
#' print(sc2)
#' 
pruneindicators<-function(x, At=0, Bt=0, sqrtIVt=0, alpha = 1.0, 
                          max.indicators=4, verbose=FALSE) {

nonnested <- function (x, selection=NULL, verbose=FALSE) {
	if(is.null(selection)) selection = rep(TRUE, nrow(x$C))
	c = x$C[selection,]
	xc = x$XC[, selection]
	combs = row.names(c)
  	keep = rep(TRUE, ncol(xc))
  	for(c1 in 1:ncol(xc)) {
  		ccx1 = xc[,c1]>0
  	  	for(c2 in 1:ncol(xc)) { 	  	
	  		ccx2 = xc[,c2]>0
  			if(c1!=c2 && keep[c2]) {
  				if(sum(ccx1 & ccx2)==sum(ccx2) && sum(ccx1)>sum(ccx2)) {
  					keep[c2] = FALSE
  					if(verbose) cat(paste(combs[c2],"nested in",combs[c1],"\n"))
  				}
  				else if(sum(ccx1 & ccx2)==sum(ccx2) && sum(ccx1 & ccx2)==sum(ccx1)) {
  					if(verbose) cat(paste(combs[c2],"equal to",combs[c1],"\n"))
  					if(sum(c[c2,])> sum(c[c1,])) keep[c2] = FALSE
  					else if(sum(c[c1,])> sum(c[c2,])) keep[c1] = FALSE
  				}
  			}
      	}
 	}	
  	return(combs[keep])
}

	initCoverage<-coverage(x)
	if(verbose) cat(paste("Coverage of initial set of ",nrow(x$C)," indicators: ", round(initCoverage*100, digits=1),"%\n", sep=""))
	
	if(length(dim(x$A))==2) {
		selection<- x$A[,"lowerCI"]>=At & x$B[,"lowerCI"]>=Bt & x$sqrtIV[,"lowerCI"]>=sqrtIVt & x$p.value <= alpha
	} else {
		selection<- x$A>=At & x$B>=Bt & x$sqrtIV>=sqrtIVt & x$p.value <= alpha
	}
    if(sum(selection)==0) {
    	if(verbose) cat(paste("No indicator is valid using the given thresholds."))
    	return()
    }
	validCoverage<-coverage(x, selection=selection)
	if(verbose) cat(paste("Coverage of valid set of ",sum(selection)," indicators: ", round(validCoverage*100, digits=1),"%\n", sep=""))
	
	if(sum(selection)>1) {
	  NN <-nonnested(x, selection=selection, verbose=FALSE)
		selection <- rownames(x$C) %in% NN
		nnCoverage<-coverage(x, selection=selection)
		if(verbose) cat(paste("Coverage of valid set of ",sum(selection)," nonnested indicators: ", round(nnCoverage*100, digits=1),"%\n", sep=""))
	
	
		c = x$C[selection,, drop = FALSE]
		group.vec = x$group.vec
		xc = x$XC[, selection]

		#Preliminaries	
  	spnames = colnames(c)
  	spplist = colnames(c)

		indnames<-rownames(c)
		k <- length(indnames)
		j=1
		continue = (!is.null(max.indicators))
    	selmodFinal = selection
   		while(continue) {
      		co <- combn(k,j) #Generate subsets of indicators
   			if(verbose) cat(paste("Checking ",ncol(co)," subsets of ", j," indicator(s)", sep=""))
	      	keep2 = rep(FALSE,ncol(co))
   	   		maxcov = 0
	      	for(coi in 1:ncol(co)) { #Check coverage of subsets of indicators
	      		if(ncol(co)>10) if(coi%%round(ncol(co)/10)==0 && verbose) cat(".")
	      		selmod = selection
	      		selmod[selection]=FALSE
	      		selmod[selection][co[,coi]]=TRUE
	      		coicov = coverage(x, selection=selmod)
	      		if(coicov>maxcov) {
	      		   bestAtPoint= coi
	      		   maxcov = coicov
	      		}
		  		keep2[coi]=coicov== nnCoverage
	      	}
	        if(verbose) cat(paste(" maximum coverage: ", round(maxcov*100, digits=1),"%\n",sep=""))
	      	if(sum(keep2)>0) { #If at least one subset has the appropriate coverage keep it
	      		best = which(keep2)[1]
	      		selmodFinal[selection]=FALSE
	      		selmodFinal[selection][co[,best]]=TRUE
	      		finalCoverage<-coverage(x, selection=selmodFinal)
				if(verbose) cat(paste("Coverage of final set of ", j, " indicators: ",round(finalCoverage*100,digits=1),"%\n", sep=""))
	      		continue = FALSE
	      	} else {
		      	if(j==max.indicators) {
	    	  		selmodFinal[selection]=FALSE
		      		selmodFinal[selection][co[,bestAtPoint]]=TRUE
					if(verbose) cat(paste("\nCoverage maximum allowed set of ", j, " indicators: ",round(maxcov*100,digits=1),"%\n", sep=""))
		      	}
	      	}
	      	if(j<min(max.indicators,k)) j= j+1
  	    	else continue = FALSE
	    }
	} else {
		if(verbose) cat(paste("One valid indicator only. Stopping.\n", sep=""))
		selmodFinal = selection
	}
    indicators2 = x
    indicators2$C = as.data.frame(subset(indicators2$C, subset=selmodFinal))
    indicators2$XC = indicators2$XC[,selmodFinal, drop=FALSE]
    if(length(dim(indicators2$A))==2) {
    	indicators2$A = indicators2$A[selmodFinal,, drop=FALSE]
	    indicators2$B = indicators2$B[selmodFinal,, drop=FALSE]
    	indicators2$sqrtIV = indicators2$sqrtIV[selmodFinal,, drop=FALSE]
    } else {
    	indicators2$A = indicators2$A[selmodFinal]
	    indicators2$B = indicators2$B[selmodFinal]
    	indicators2$sqrtIV = indicators2$sqrtIV[selmodFinal]
    }
    indicators2$p.value = indicators2$p.value[selmodFinal]

	return(indicators2)
}



