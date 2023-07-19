
#' Resource niche metrics for a single niche
#' 
#' Function \code{nichepref} computes the species resource preference from a the 
#' species resource use (and resource availability when given). Function \code{nichecentroid} 
#' computes the centroid on the resource space for a set of species. Function \code{nichevar} 
#' computes the multivariate resource variance for a set of species (i.e. niche breadth). 
#' In all functions resources are given in distance matrix \code{D}, the species resource use is given in \code{P} and the availability of resources, if present, are given by vector \code{q}. 
#' 
#' @param P Data frame containing the relative or absolute usage that a set of species (in rows) make of a set of resources (in columns)
#' @param D Object of type \code{\link{dist}} containing distance values between resources. If no distance matrix is provided (i.e. if \code{D==NULL}), the distances between resources is assumed to be maximum.
#' @param q Vector with the availability of each resource
#' @param mode Either \code{mode = "single"} (rows of matrix P are individual observations to be pooled for a single niche) or \code{mode = "multiple"} (rows in P represent different niches)
#' @param Np Vector with the number of observations per species from which the values in \code{P} come (in \code{mode = "multiple"})
#' @param Nq The number of observations per species from which the values in \code{q} come (in \code{mode = "multiple"})
#' @param nboot Number of boostrap samples used to compute bias-corrected percentile confidence intervals
#' @param alpha Used to set the confidence level (i.e. \code{alpha = 0.05} means 95 percent confidence interval)
#'
#' @details
#' The method is described in De Caceres et al. (2010). If the distance matrix is not specified 
#' (i.e. if \code{D=NULL}) the function assumes that all resources are at a maximum distance (\code{d=1}).
#'  If the resource availability vector \code{q} is given then the values in \code{P} are taken as 
#'  assessments of resource use and the species preference is calculated taking into account resource 
#'  availability. Otherwise resource use is equated to resource preference. Moreover, most functions can 
#'  compute bootstrap confidence intervals following the bias-corrected percentile method (Manly 2007). 
#'  If \code{mode = "multiple"} and \code{Np != NULL}, bootstrap samples for a given species are generated
#'  assuming a multinomial distribution with the proportions calculated from the corresponding row values 
#'  in \code{P}, and the number of observations comes from the corresponding element in \code{Np}. 
#'  If \code{mode = "single"} then the bootstrapped units are the rows of matrix \code{P}. In both cases, 
#'  if \code{Nq} is indicated the availability of resources is also bootstrapped. 
#'  The bias-corrected percentile method was described for overlap niche measures 
#'  in Mueller and Altenberg (1985) and is extended here for all niche metrics.
#' 
#' @references 
#' Mueller, L.D. and L. Altenberg. 1985. Statistical Inference on Measures of Niche Overlap. Ecology 66:1204-1210.
#' 
#' Manly, B.F.J. 2007. Randomization, bootstrap and Monte Carlo methods in biology. Chapman and Hall texts in statistical science series. 2nd edition.
#' 
#' De Caceres, M., Sol, D., Lapiedra, O. and P. Legendre. (2011) A framework for estimating niche metrics using the resemblance between qualitative resources. Oikos 120: 1341-1350.
#' 
#' @author Miquel De Caceres Ainsa, EMF-CREAF
#' 
#' @seealso See \code{\link{nicheoverlap}} for descriptors comparing two niches.
#' 
#' @return
#' Function \code{nichepref} returns a matrix of species relative preference. 
#' Function \code{nichevar} returns a vector with the variance of the resources 
#' used for each species in \code{P}. Function \code{nichecentroid} returns a matrix niche centroid in the resource space for each species in \code{df}. If bootstrap confidence intervals are asked then the three functions also compute two extra data containing respectively the lower and upper bounds of the confidence intervals obtained following the bias-corrected percentile method. Function \code{nichearea} returns the area of the convex hull occupied by the resources used for each species in \code{P}. 
#' 
#' @export
#'
#' @name nichevar
#' @examples
#' # Loads example data
#' data(birds)
#' 
#' # The niche metrics using distances among resources and assuming equal availability of resources
#' nichepref(birdsbreed, D = resourceD) 
#' nichevar(birdsbreed, D = resourceD) 
#' nichecentroid(birdsbreed, D = resourceD) 
#' 
#' # The niche metrics using distances among resources and computes 
#' # 95 percent confidence intervals
#' nichepref(birdsbreed, D = resourceD, mode="multiple", 
#'           Np = rowSums(birdsbreed), Nq = 100) 
#' nichevar(birdsbreed, D = resourceD, mode="multiple", 
#'          Np = rowSums(birdsbreed), Nq = 100) 
#' nichecentroid(birdsbreed, D = resourceD, mode="multiple", 
#'               Np = rowSums(birdsbreed), Nq = 100) 
#' 
#' # Same computations with different resource availability
#' nichepref(birdsbreed, D = resourceD, 
#'           q = c(0.18, 0.24, 0.22, 0.21, 0.15), mode="multiple")
#' nichevar(birdsbreed, D = resourceD, 
#'          q = c(0.18, 0.24, 0.22, 0.21, 0.15), mode="multiple")
#' nichecentroid(birdsbreed, D = resourceD, 
#'               q = c(0.18, 0.24, 0.22, 0.21, 0.15), mode="multiple")
#' 
#' # The niche metrics using distances among resources and 
#' # computes 95 percent confidence intervals
#' nichepref(birdsbreed, D = resourceD, 
#'           q = c(0.18, 0.24, 0.22, 0.21, 0.15), mode="multiple", Np = rowSums(birdsbreed), Nq = 100)
#' nichevar(birdsbreed, D = resourceD, 
#'          q = c(0.18, 0.24, 0.22, 0.21, 0.15), mode="multiple", Np = rowSums(birdsbreed), Nq = 100)
#' nichecentroid(birdsbreed, D = resourceD, 
#'               q = c(0.18, 0.24, 0.22, 0.21, 0.15), mode="multiple",  Np = rowSums(birdsbreed), Nq = 100)
#' 
nichevar <- function (P, D = NULL, q = NULL, mode="multiple", Np = NULL, Nq = NULL, nboot = 1000, alpha=0.05) {
    if (!inherits(P, "data.frame")) stop("Non convenient dataframe for species resource use")
    if (!is.null(D)) {
        if (!inherits(D, "dist")) stop("Object of class 'dist' expected for distance")
        D <- as.matrix(D)
        if (ncol(P) != nrow(D)) stop("The number of columns in P must be equal to the number of items in D")
        D <- as.dist(D)
    }
    if(!is.null(Np) && mode=="multiple") {
    	   if(length(Np)!=nrow(P)) stop("The number of items in Np must be equal to the number of rows in P")
    }
    if(!is.null(q)) {
    	   if(length(q)!=ncol(P)) stop("The number of items in q must be equal to the number of columns in P")
    	   q = q/sum(q) #Just to check that we have proportions
    } else {
        #If no availability is provided, then all resources are equally available
    	   q = rep(1/ncol(P),ncol(P))
    }
    	
    #If no distance matrix is provided, the distance between resources is assumed to be maximum
    if (is.null(D)) D <- as.dist((matrix(1, ncol(P), ncol(P)) - diag(rep(1, ncol(P)))))
    
    # Computes the niche breadth from the resource preferences of the target and the resource relationships
    nichevar1<-function(f, D) {
       if (is.na(sum(f))) v <- NA
       else if (sum(f) < 1e-16) v <- 0
       else v <- (f %*% (as.matrix(D)^2) %*% f)/(2*(sum(f)^2))
   		return(v)
    	}
    	
    	# Returns preference from a resource use vector (considering resource availability in desired)
    	getF<-function(p,q=NULL) {
    		if(!is.null(q)) {
    			a = p/q
    			return(a/sum(a))
    		} else { #Just to check that we have proportions
    			return(p/sum(p))
    		}
    	}
    	
    if(!is.null(Np) || mode=="single") nc = 3
    else nc = 1
    
    #Rows in P are different niches
    if(mode=="multiple") {
	    B <- as.data.frame(matrix(0,nrow=nrow(P), ncol=nc))
 	   rownames(B) <- row.names(P)
 	   for (i in 1:nrow(P)) {
 	   	  pi = as.numeric(P[i,])
 	   	  B[i,1] = nichevar1(getF(pi,q), D)
 	   	  
 	   	  if(!is.null(Np)) {
  	  		  		BB = vector("numeric",length=nboot)
 	   		  	if(sum(is.na(getF(pi)))==0) {
	    		  	   #Generate bootstrap samples from multinomial distribution
 	   		  	   psamp = rmultinom(nboot,Np[i],getF(pi))
 	   		  	   if(!is.null(Nq)) qsamp = rmultinom(nboot,Nq,q)
 	   		  	   for(b in 1:nboot) {
 	   		  	   	 if(!is.null(Nq)) BB[b] = nichevar1(getF(psamp[,b],qsamp[,b]),D)
 	   		  	   	 else BB[b] = nichevar1(getF(psamp[,b],q),D)
 	   		  	   }
  	  		  	   	#Some NA may appear because of zeroes in qsamp
	    		  	   BB = BB[!is.na(BB)]
				    	 #Compute Bias-corrected percentile method (Manly 2007: pp52-56)
 	   		  	   z0 = qnorm(sum(BB<B[i,1])/length(BB))
 	   		  	   lj = floor(length(BB)*pnorm(2*z0+qnorm(alpha/2)))
 	   		  	   uj = floor(length(BB)*pnorm(2*z0+qnorm(1-(alpha/2))))
 	   		  	   if(lj > 0 && uj > 0) {
		 	   	  	   sbb = sort(BB)
		 	   	  	   B[i,2] = sbb[lj]
			   	  	   B[i,3] = sbb[uj]
		    	  	   }
 	   	  	  }
 	   	  	}
 	   }
 	  }
    #Rows in P are observations
    else if(mode=="single") {
	   B <- as.data.frame(matrix(0,nrow=1, ncol=nc))
 	   rownames(B) <- "Niche"
 	   B[1,1] = nichevar1(getF(colSums(P),q), D)
 	   	  
  	  	 BB = vector("numeric",length=nboot)
 	   if(!is.null(Nq)) qsamp = rmultinom(nboot,Nq,q)
  	  	 for(b in 1:nboot) {
	 	   psamp = colSums(P[sample(1:nrow(P),replace=TRUE),])
	 	   if(!is.null(Nq)) {
 		   	 BB[b] = nichevar1(getF(psamp,qsamp[b]),D)
 		   } else {
 		   	 BB[b] = nichevar1(getF(psamp,q),D)
 		   }
  	  	 }
  	  	#Some NA may appear because of zeroes in qsamp
	   BB = BB[!is.na(BB)]
	   #Compute Bias-corrected percentile method (Manly 2007: pp52-56)
 	   z0 = qnorm(sum(BB<B[1,1])/length(BB))
 	   lj = floor(length(BB)*pnorm(2*z0+qnorm(alpha/2)))
 	   uj = floor(length(BB)*pnorm(2*z0+qnorm(1-(alpha/2))))
 	   if(lj > 0 && uj > 0) {
		 	  sbb = sort(BB)
		 	  B[1,2] = sbb[lj]
			  B[1,3] = sbb[uj]
		 }
 	  } 
 	  if(nc==1) names(B) <- "B"
 	  else names(B) <- c("B","LC", "UC")
    return(B)
}