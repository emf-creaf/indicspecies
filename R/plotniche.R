#' Draws a single resource niche
#' 
#' Function \code{plotniche} draws a resource niche on the resource space computed by means of principal coordinates analysis. Resource relationships are given in distance matrix \code{D}, the species resource use is given in \code{P} and the availability of resources is given by vector \code{q}. 
#'
#' @param P Data frame containing the relative or absolute usage that a set of species (in rows) make of a set of resources (in columns).
#' @param D Object of type \code{\link{dist}} containing distance values between resources. If no distance matrix is provided (i.e. if \code{D==NULL}), the distances between resources is assumed to be maximum.
#' @param q Vector with the availability of each resource.
#' @param mode Either \code{mode = "single"} (rows of matrix P are individual observations to be pooled for a single niche) or \code{mode = "multiple"} (rows in P represent different niches).
#' @param Np Vector with the number of observations per species from which the values in \code{P} come (in \code{mode = "multiple"}).
#' @param Nq The number of observations per species from which the values in \code{q} come.
#' @param nboot Number of boostrap samples used to compute bias-corrected percentile confidence intervals.
#' @param alpha Used to set the confidence level (i.e. \code{alpha = 0.05} means 95 percent confidence interval).
#' @param species Specifies which species niche is to be plot. This parameter is mandatory and can be either an numeric index or a string for a species name.
#' @param axes PCoA axes used for plotting.
#' @param chull Logical flag indicating whether or not convex hulls should be drawn (only in \code{type="single"}).
#' @param bubbles Logical flag to draws bubbles proportional to resource preference data.
#' @param writeName Logical flag indicating whether or not the name of the species should be drawn beside the centroid.
#' @param add If \code{TRUE}, the current plot is used. This is helpful to draw more than one species on the same plot (see examples). 
#' @param col Color of the centroid and confidence interval arrows.
#' @param lty Line type of the confidence interval arrows.
#' @param ... Additional graphical parameters.
#'
#' @details
#' The method is described in De Caceres et al. (in prep). If the distance matrix is not specified (i.e. if \code{D=NULL}) the function assumes that all resources are at a maximum distance (\code{d=1}). If the resource availability vector \code{q} is given then the values in \code{P} are taken as assessments of resource use and the species preference is calculated taking into account resource availability. Otherwise resource use is equated to resource preference. The function can also plot bootstrap confidence intervals following the bias-corrected percentile method (Manly 2007). If If \code{mode = "multiple"} and \code{Np != NULL}, bootstrap samples for a given species are generated assuming a multinomial distribution with the proportions calculated from the corresponding row values in \code{P}, and the number of observations comes from the corresponding element in \code{Np}. If \code{mode = "single"} then the bootstrapped units are the rows of matrix \code{P}. In both cases, if \code{Nq} is indicated, the availability of resources is also bootstrapped. The bias-corrected percentile method was described for overlap niche measures in Mueller and Altenberg (1985) and is extended here for all niche metrics except \code{nichearea}.
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
#' @seealso See \code{\link{nichevar}} and \code{\link{nicheoverlap}} to obtain a resource niche metrics.
#' 
#' @export
#'
#' @examples
#' # Loads example data
#' data(birds)
#'
#' plotniche(birdsbreed, D = resourceD, mode="multiple", species=10) 
#' plotniche(birdsbreed, D = resourceD, mode="multiple", 
#'           Np = rowSums(birdsbreed), Nq = 100, species=10) 
#' plotniche(birdsbreed, D = resourceD, 
#'           q = c(0.18, 0.24, 0.22, 0.21, 0.15), mode="multiple", species=10)
#' plotniche(birdsbreed, D = resourceD, 
#'           q = c(0.18, 0.24, 0.22, 0.21, 0.15), mode="multiple", 
#'           Np = rowSums(birdsbreed), Nq = 100, species=10)
#' 
#' #draw two species
#' plotniche(birdsbreed, D = resourceD, mode="multiple", 
#'           Np = rowSums(birdsbreed), Nq = 100, species=10, writeName=TRUE, 
#'           bubbles=FALSE, chull=FALSE) 
#' plotniche(birdsbreed, D = resourceD,  mode="multiple", 
#'           Np = rowSums(birdsbreed), Nq = 100, species=1, writeName=TRUE, 
#'           bubbles=FALSE, chull=FALSE, add=TRUE, col="red", lty=2) 
plotniche <-
function(P, D = NULL, q = NULL, mode="multiple", Np = NULL, Nq = NULL, nboot = 1000, alpha=0.05, species=NULL, axes=c(1,2), 
         chull=TRUE, bubbles=TRUE, writeName=FALSE, add=FALSE, col="black", lty=1,...) {
	 
    #If no distance matrix is provided, the distance between resources is assumed to be maximum
    if (is.null(D)) D <- as.dist((matrix(1, ncol(P), ncol(P)) - diag(rep(1, ncol(P)))))

	 #Computes metric MDS
	 cmd = cmdscale(D,eig=TRUE,k= ncol(P)-1, add = TRUE)
	
	 eigp = 100*cmd$eig/sum(cmd$eig)
	
	if(!add) {
		plot(cmd$points[,axes], xlab=paste("PCoA ",axes[1]," (",format(eigp[axes[1]],digits=3),"%)",sep=""), ylab=paste("PCoA ",axes[2]," (",format(eigp[axes[2]],digits=3),"%)",sep=""), cex=0.2,...)
		text(cmd$points[,axes],labels=names(P), cex=1.0, pos=3, offset=0.3)
	}
	if(mode=="multiple") {
	 	if(is.null(species)) stop("Please, provide a species to be plot")
	 	if(is.numeric(species)) {
				isp = species
				spname = row.names(P)[isp]
	 	 }	
	 	 else if(is.character(species)) {
			 	spname = species	
				isp = which(row.names(P)==spname)
		 }
    	 # Preference from a resource use vector (considering resource availability in desired)
		 if(isp<0) stop("Species not found")
		 cat(paste("\n Plotting species :", spname, " #:",isp),"\n\n")
	}	
		
	centr = nichecentroid(P=P, D=D, q = q, mode=mode, Np = Np, Nq = Nq, nboot = nboot, alpha=alpha)
	pref = nichepref(P=P, D=D, q = q, mode=mode, Np = Np, Nq = Nq, nboot = nboot, alpha=alpha)
		
	#Bubbles proportional to resource preference
	if(bubbles) {
			if(mode=="multiple") { 
				if(!is.null(Np)) {
					a = subset(cmd$points[,axes], pref$F[isp,]>0)
					symbols(a, circles=(pref$F[isp,pref$F[isp,]>0]/20), inches=FALSE,add=TRUE,lwd=2)
					symbols(a, circles=(pref$LF[isp,pref$F[isp,]>0]/20), inches=FALSE,add=TRUE,lwd=1, lty=2)
					symbols(a, circles=(pref$UF[isp,pref$F[isp,]>0]/20), inches=FALSE,add=TRUE,lwd=1, lty=2)
				} else {
					a = subset(cmd$points[,axes], pref[isp,]>0)
					symbols(a, circles=(pref[isp,pref[isp,]>0]/20), inches=FALSE,add=TRUE,lwd=2)
				}
			} else if(mode=="single") {
				a = subset(cmd$points[,axes], pref[1,]>0)
				symbols(a, circles=(pref[1,pref[1,]>0]/20), inches=FALSE,add=TRUE,lwd=2)
				symbols(a, circles=(pref[2,pref[1,]>0]/20), inches=FALSE,add=TRUE,lwd=1, lty=2)
				symbols(a, circles=(pref[3,pref[1,]>0]/20), inches=FALSE,add=TRUE,lwd=1, lty=2)
			}
		}
		
		if(mode=="multiple") {
			#Resource niche centroid (with confidence intervals)
		  if(!is.null(Np)) {
				points(centr$C[isp,axes],pch=21, bg = col, col=col)
				arrows(x0=centr$LC[isp,axes[1]],y0=centr$C[isp,axes[2]],x1=centr$UC[isp,axes[1]], angle=90, code=3, length=0.1, col=col, lty=lty)
				arrows(y0=centr$LC[isp,axes[2]],x0=centr$C[isp,axes[1]],y1=centr$UC[isp,axes[2]], angle=90, code=3, length=0.1, col=col, lty=lty)
				if(writeName) text(x=centr$C[isp,axes[1]], y=centr$C[isp,axes[2]], pos=3, labels = spname, col=col)
			}
			else {
				points(centr[isp,axes],pch=21, bg = col, col=col)
				if(writeName) text(x=centr[isp,axes[1]], y=centr[isp,axes[2]], pos=3, labels = spname, col=col)
			}
		}
		else if(mode=="single") {
				points(centr[1, axes],pch=21, bg = col, col=col)
				arrows(x0=centr[2,axes[1]],y0=centr[1,axes[2]],x1=centr[3,axes[1]], angle=90, code=3, length=0.1, col=col, lty=lty)
				arrows(y0=centr[2,axes[2]],x0=centr[1,axes[1]],y1=centr[3,axes[2]], angle=90, code=3, length=0.1, col=col, lty=lty)
				if(writeName) text(x=centr[1,axes[1]], y=centr[1,axes[2]], pos=3, labels = species, col=col)
			
		}
		
		#Convex hull
		if(chull) {
			if(mode=="multiple") {
				if(!is.null(Np)) a = subset(cmd$points[,axes], pref$F[isp,]>0)
				else a = subset(cmd$points[,axes], pref[isp,]>0)
			}
			else if(mode=="single") a = subset(cmd$points[,axes], pref[1,]>0)
		  chp = chull(a[,1],a[,2])
			polygon(a[chp,1],a[chp,2],lwd=1)
		}
		
}
