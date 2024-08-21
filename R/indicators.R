#' Indicator analysis for species combinations
#'
#' Determines the indicator value of species combinations.
#' 
#' @param X A community data table with sites in rows and species in columns. This table can contain either presence-absence or abundance data.
#' @param cluster A vector containing the classification of sites into site groups.
#' @param group The label corresponding to the target site group. 
#' @param func The indicator value variant to be used, either "IndVal" (non-equalized) or "IndVal.g" (group-equalized).
#' @param min.order Minimum number of species conforming species combinations.
#' @param max.order Maximum number of species conforming species combinations.
#' @param max.indicators Maximum number of valid indicators to be kept. If \code{NULL}, then all valid indicators are kept.
#' @param At Threshold for positive predictive value used to select valid indicators. Combinations with lower values are not kept.
#' @param Bt Threshold for sensitivity used to select valid indicators. Combinations with lower values are not kept.
#' @param sqrtIVt Threshold for (square root of) indicator value. Combinations with lower values are not kept.
#' @param alpha Threshold for statistical significance of indicator value. Combinations with higher p-values are not kept.
#' @param control A list of control values describing properties of the permutation test design, as returned by a call to \code{\link[permute]{how}}.
#' @param permutations A custom matrix of permutations, to be used if \code{control = NULL}, with permutations in rows and site indices in columns.
#' @param print.perm If TRUE, prints permutation numbers after each set of 100 permutations.
#' @param nboot.ci Number of bootstrap samples for confidence intervals. If \code{nboot.ci = NULL} then confidence intervals are not estimated.
#' @param alpha.ci Error in confidence intervals.
#' @param XC If TRUE, outputs the abundance/occurrence matrix of species combinations.
#' @param enableFixed If TRUE, uses species that occur in all sites as fixed elements and creates combinations with the remaining ones.
#' @param verbose If TRUE, prints the results of each step.
#'
#' @details
#' Function \code{indicators} creates explores the indicator value of the simultaneous occurrence of sets of species (i.e. species combinations). The method is described in De \enc{Cáceres}{Caceres} et al. (2012) and is a generalization of the Indicator Value method of \enc{Dufrêne}{Dufrene} & Legendre (1997). The minimum and maximum number of species conforming the species combination can be controlled using \code{min.order} or \code{max.order}. For each combination of species it determines its positive predictive value (A), sensitivity (B) and the square root of indicator value (sqrtIV). Statistical significance of indicators for the target site group is determined by internal calls to function \code{\link{signassoc}}. Additionally, if \code{nboot.ci} is not null then bootstrap confidence intervals are determined with the specified \code{alpha} level, as explained in De \enc{Cáceres}{Caceres} & Legendre (2009). The combinations to be kept can be restricted to those whose positive predictive value, sensitivity and/or indicator value are equal or greater than input thresholds. Function \code{print} allows printing the results in a nice table, whereas \code{summary} provides information about candidate species, combinations and coverage of the set of indicators. Function \code{plot} draws the statistics against the order (i.e. the number of species) of the combination.
#' 
#' @return
#' An object of class \code{indicators} with:
#' \item{candidates}{The vector of initial candidate species.}
#' \item{finalsplist}{The vector of species finally selected for combinations.}
#' \item{C}{A matrix describing all the combinations studied.}
#' \item{XC}{A matrix containing the abundance/occurrence of each species combination.}
#' \item{A}{Positive predictive power of species combinations. If \code{nboot} is not missing then this includes the lower and upper bounds of the confidence interval.}
#' \item{B}{Sensitivity of species combinations. If \code{nboot} is not missing then this includes the lower and upper bounds of the confidence interval.}
#' \item{sqrtIV}{Square root of indicator value of species combinations. If \code{nboot} is not missing then this includes the lower and upper bounds of the confidence interval.}
#' \item{sign}{P-value of the permutation test of statistical significance.}
#' \item{group.vec}{A logical vector indicating the membership to the target group.}
#' 
#' @references 
#' De \enc{Cáceres}{Caceres}, M., Legendre, P., Wiser, S.K. and Brotons, L. 2012. Using species combinations in indicator analyses. Methods in Ecology and Evolution 3(6): 973-982.
#' 
#' De \enc{Cáceres}{Caceres}, M. and Legendre, P. 2009. Associations between species and groups of sites: indices and statistical inference. Ecology 90(12): 3566-3574.
#' 
#' \enc{Dufrêne}{Dufrene}, M. and P. Legendre. 1997. Species assemblages and indicator species: The need for a flexible asymetrical approach. Ecological Monographs 67:345-366.
#' 
#' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, EMF-CREAF
#' 
#' @seealso \code{\link{predict.indicators}},\code{\link{pruneindicators}}, \code{\link{coverage}}, \code{\link{multipatt}}, \code{\link{strassoc}}, \code{\link{signassoc}}
#' 
#' @export
#'
#' @name indicators
#' 
#' @examples
#' library(stats)
#' 
#' data(wetland) ## Loads species data
#' 
#' ## Creates three clusters using kmeans
#' wetkm <- kmeans(wetland, centers=3) 
#' 
#' ## Number of sites in each group
#' table(wetkm$cluster)
#' 
#' ## Run indicator analysis with species combinations for the first group
#' sc <- indicators(X=wetland, cluster=wetkm$cluster, group=1, verbose=TRUE, 
#'                  At=0.5, Bt=0.2)
#' 
#' #Prints the results
#' print(sc)
#' 
#' ## Plots positive predictive power and sensitivity against the order of 
#' ## combinations
#' plot(sc, type="A")
#' plot(sc, type="B")
#' 
#' ## Run indicator analysis with species combinations for the first group, 
#' ## but forcing 'Orysp' to be in all combinations
#' sc2 <- indicators(X=wetland, cluster=wetkm$cluster, group=1, verbose=TRUE, 
#'                   At=0.5, Bt=0.2, enableFixed=TRUE)
#'                   
indicators <- function (X, cluster, group, func="IndVal", min.order = 1, max.order = 5, max.indicators=NULL, 
                        At = 0, Bt=0, sqrtIVt =0, 
                        control = how(), permutations = NULL, print.perm = FALSE,
                        nboot.ci=NULL, alpha.ci=0.05, XC = TRUE, enableFixed = FALSE, verbose=FALSE) {
	                 
  # Turn into a matrix (if not)
  X = as.matrix(X)
  
  func <- match.arg(func, c("IndVal", "IndVal.g"))                                                                                                             
  if(sum(is.na(cluster))>0) stop("Cannot deal with NA values. Remove and run again.")
  if(sum(is.na(X))>0) stop("Cannot deal with NA values. Remove and run again.")
  if(min.order<1) stop("min.order must be at least 1.")
  if(min.order>max.order) stop("min.order must be equal or larger than max.order.")
  nsites = nrow(X)

  cluster= as.factor(cluster)
  group <- as.character(group)
  group <- match.arg(group, levels(cluster))                                                                                                             
  if(verbose) cat(paste("Target site group: ",group,"\n", sep=""))
  group.vec <- cluster==group
  group.vec[is.na(group.vec)] = FALSE

  #Get species names
  spplist = colnames(X)
  if(verbose) cat(paste("Number of candidate species: ",length(spplist),"\n", sep=""))
  if(length(spplist)==1) stop("At least two species are necessary.")
  
  #Select rows that contain the species or the group
  ng = sum(group.vec)
  if(verbose) cat(paste("Number of sites:",nsites,"\n"))
  if(verbose) cat(paste("Size of the site group:",ng,"\n"))
  
  #Study frequency
  freq = colSums(ifelse(X[group.vec,]>0,1,0))/sum(group.vec)
  if(enableFixed) {
    fixedSpecies = (freq==1)
  } else {
    fixedSpecies = rep(FALSE, length(spplist))
  }
  numFixed = sum(fixedSpecies)
  if(verbose & enableFixed) cat(paste("Number of fixed species:",numFixed,"\n"))
  fixedPos = which(fixedSpecies)
  if(verbose & numFixed>0) print(fixedPos)
  

  evalCombs<-function(spvec, dvec, verbose=FALSE) {
    comblist<-vector("list",0)
    sc.ab<-apply(X[,spvec, drop=FALSE],1,min)
    if(sum(sc.ab)>0) {
      scg = sc.ab[group.vec]
      if(func=="IndVal.g") {
        mg = (sum(scg)/ng)
        A = mg/sum(tapply(sc.ab,cluster, "mean"))
      } else {
        A = sum(scg)/sum(sc.ab)
      }
      B = sum(scg>0)/ng
      sqrtIV = sqrt(A*B)
      if(A>=At & B>=Bt & sqrtIV>=sqrtIVt) {  
        if(length(spvec)>=min.order) comblist<-c(comblist,list(spvec))
      }
      if(B>Bt && (length(spvec)-numFixed)<max.order && length(dvec)>0) {##If B is not too small and order is allowed we can explore further
        for(j in 1:length(dvec)) {
          if(verbose) cat(paste("Starting species ",dvec[j],"..."))          
          comblist<-c(comblist,evalCombs(c(spvec,dvec[j]), dvec[-(1:j)], verbose=FALSE))
          if(verbose) cat(paste(" accepted combinations:",length(comblist),"\n"))
        }      
      }
    } 
    return(comblist)
  }
  k = length(spplist)-numFixed #Number of species to combine
  if(verbose & enableFixed) cat(paste("Number of species to combine: ",k,"\n", sep=""))
  veck = (1:length(spplist))[!fixedSpecies] #Vector of species indices to combine
  
  comblistDef<-vector("list",0)
  if(length(fixedPos)>0) {
    comblistDef<-c(comblistDef,evalCombs(fixedPos, veck,verbose=TRUE))
  } else {
    for(i in 1:length(veck)) {
      if(verbose) cat(paste("Starting species ",veck[i],"..."))
      comblistDef<-c(comblistDef,evalCombs(c(veck[i],fixedPos), veck[-(1:i)], verbose=FALSE))
      if(verbose) cat(paste(" accepted combinations:",length(comblistDef),"\n"))
    }    
  }
  
  #Create structures to store data
  nc = length(comblistDef)
  if(verbose) cat(paste("Number of valid combinations: ",nc,"\n", sep=""))
  if(nc==0) return()
  trim =FALSE
  if(!is.null(max.indicators)) {
    if(nc>max.indicators) {
      if(verbose) cat(paste("Maximum number of valid combinations exceeded.\n", sep=""))    
      trim = TRUE
    }
  }
  Astat = numeric(nc)
  Bstat = numeric(nc)
  if(trim) {
    for(i in 1:nc) {
      spvec = as.numeric(comblistDef[[i]])
      sc.ab <-apply(X[,spvec, drop=FALSE],1,min)
      scg = sc.ab[group.vec]
      if(func=="IndVal.g") {
        mg = (sum(scg)/ng)
        Astat[i] = mg/sum(tapply(sc.ab,cluster, "mean"))
      } else {
        Astat[i] = sum(scg)/sum(sc.ab)
      }
      Bstat[i] = sum(scg>0)/ng
    }
    sqrtIVstat = sqrt(Astat*Bstat)    
    sel = order(sqrtIVstat,decreasing = TRUE)[1:max.indicators]
    Astat = Astat[sel]
    Bstat = Bstat[sel]
    sqrtIVstat = sqrtIVstat[sel]
    nc = max.indicators
    Cvalid<-matrix(0,nrow=nc,ncol=length(spplist))
    colnames(Cvalid)<-spplist
    XC = matrix(0, nrow=nrow(X), ncol=nc)
    for(i in 1:nc) {
      spvec = as.numeric(comblistDef[[sel[i]]])
      XC[,i] <-apply(X[,spvec, drop=FALSE],1,min)
      Cvalid[i,spvec] = 1
    }
  } else {
    Cvalid = matrix(0,nrow=nc,ncol=length(spplist))
    XC = matrix(0, nrow=nrow(X), ncol=nc)
    colnames(Cvalid)<-spplist
    for(i in 1:nc) {
      spvec = as.numeric(comblistDef[[i]])
      sc.ab <-apply(X[,spvec, drop=FALSE],1,min)
      XC[,i]<-sc.ab
      scg = sc.ab[group.vec]
      if(func=="IndVal.g") {
        mg = (sum(scg)/ng)
        Astat[i] = mg/sum(tapply(sc.ab,cluster, "mean"))
      } else {
        Astat[i] = sum(scg)/sum(sc.ab)
      }
      Bstat[i] = sum(scg>0)/ng
      Cvalid[i,spvec] = 1
    }
    sqrtIVstat = sqrt(Astat*Bstat)    
  }
  cn = 1:ncol(XC)
  for(i in 1:ncol(XC)) cn[i] = paste(spplist[Cvalid[i,]==1], collapse="+")
  rownames(Cvalid) = cn
  rownames(XC) = rownames(X)
  
  #Remove species that do not appear in any valid combination
  selSpp = colSums(Cvalid)>0
  if(verbose) cat(paste("Number of remaining species:",sum(selSpp),"\n"))
  Cvalid = Cvalid[,selSpp]
  nspp <- sum(selSpp)
  
  # Calculate statistical significance by calling signassoc
  if(verbose) cat(paste("Calculating statistical significance (permutational test)...\n"))
  mode = ifelse(func=="IndVal.g",1,0)
  p.value = signassoc(XC, cluster = cluster, mode = mode, control = control, 
                      permutations = permutations, print.perm = print.perm)[,group]
  
  # Calculate bootstrap confidence intervals for sensitivity and ppp of valid combinations
  if(!is.null(nboot.ci)) {
  	  if(nboot.ci<100) nboot.ci=100 #Minimum of 100 bootstrap replicates
	  if(verbose) {
  			cat(paste("Calculating bootstrap confidence intervals"))
  	  }
	  dmbA = matrix(NA,nrow=nboot.ci,ncol=nc)
	  dmbB = matrix(NA,nrow=nboot.ci,ncol=nc)
	  dmbIV = matrix(NA,nrow=nboot.ci,ncol=nc)
	  for(b in 1:nboot.ci) {
	  	  if(b%%round(nboot.ci/10)==0 && verbose) cat(".")
		  bi = sample(nsites,replace=TRUE)
		  ngb = sum(group.vec[bi])
		  XCB = as.matrix(XC[bi,])
		  XCBg = as.matrix(XCB[group.vec[bi],])
	  	  if(func=="IndVal.g") {
	  	  	kk <- colSums(apply(XCB,MARGIN=2,FUN=tapply,cluster[bi],"mean"))
		  	dmbA[b,] = colMeans(XCBg)/kk
  		  } else {
  			dmbA[b,] = colSums(XCBg)/colSums(XCB)
  		  }
		  dmbB[b,] = colSums(as.matrix(ifelse(XCBg>0,1,0)))/ngb
		  dmbIV[b,] = sqrt(dmbA[b,]*dmbB[b,])
	  }
	  if(verbose) cat(paste("\n"))
	  dmlowerA = rep(0,nc)
	  dmupperA = rep(0,nc)
	  dmlowerB = rep(0,nc)
	  dmupperB = rep(0,nc)
	  dmlowerIV = rep(0,nc)
	  dmupperIV = rep(0,nc)
	  for(i in 1:nc) {	
			sdmb = sort(dmbA[,i])			
			dmlowerA[i]=sdmb[(alpha.ci/2.0)*nboot.ci]
			dmupperA[i]=sdmb[(1-(alpha.ci/2.0))*nboot.ci]
			sdmb = sort(dmbB[,i])
			dmlowerB[i]=sdmb[(alpha.ci/2.0)*nboot.ci]
			dmupperB[i]=sdmb[(1-(alpha.ci/2.0))*nboot.ci]
			sdmb = sort(dmbIV[,i])
			dmlowerIV[i]= sdmb[(alpha.ci/2.0)*nboot.ci]
			dmupperIV[i]= sdmb[(1-(alpha.ci/2.0))*nboot.ci]
	  }
	  sA = cbind(Astat,dmlowerA,dmupperA)
  	  colnames(sA) = c("stat", "lowerCI", "upperCI")
  	  rownames(sA) = row.names(Cvalid)
	  sB = cbind(Bstat,dmlowerB,dmupperB)
  	  colnames(sB) = c("stat", "lowerCI", "upperCI")
  	  rownames(sB) = row.names(Cvalid)
	  sIV = as.data.frame(cbind(sqrtIVstat,dmlowerIV,dmupperIV))
  	  names(sIV) = c("stat", "lowerCI", "upperCI")
  	  rownames(sIV) = row.names(Cvalid)
  
  } else {
  	  sA = Astat
  	  sB = Bstat
  	  sIV = sqrtIVstat
  }
  
  result = list(func = func, X = X, cluster = cluster, group.vec =group.vec, candidates = spplist, finalsplist= spplist[selSpp], 
                C=Cvalid, XC=XC, A=sA, B=sB, sqrtIV=sIV, p.value = p.value)
  class(result) = "indicators"
  return(result)
}


#' @rdname indicators
#' @param x An object of class 'indicators'.
#' @param selection A logical vector used to restrict, a priori, the species combinations to be printed.
#' @param confint Flag to indicate that confidence interval bounds are desired when printing.
print.indicators <- function (x, At = 0, Bt = 0, sqrtIVt = 0, alpha = 1.0, 
                              selection = NULL, confint=FALSE,...) {
  if(is.null(selection)) selection = rep(TRUE, nrow(x$C))
  if(length(dim(x$A))==2) {
    A = as.data.frame(x$A[selection, , drop = FALSE])
    B = as.data.frame(x$B[selection, , drop = FALSE])
    sqrtIV = as.data.frame(x$sqrtIV[selection, , drop = FALSE])
  } else {
    A = x$A[selection]
    B = x$B[selection]
    sqrtIV = x$sqrtIV[selection]
  }
  p.value = x$p.value[selection]
  C = subset(x$C,selection)
  spnames = colnames(C)	
  if(length(dim(A))==2) {
    sel = A$stat>=At & B$stat>=Bt & sqrtIV$stat>=sqrtIVt & p.value<=alpha
    if(confint) {
      nc= 10
      A = A[sel,]
      B = B[sel,]
      sqrtIV = sqrtIV[sel,]
    } else {
      nc= 4	
      A = A$stat[sel]
      B = B$stat[sel]
      sqrtIV = sqrtIV$stat[sel]
    }
    p.value = p.value[sel]
  } else {
    sel = A>=At & B>=Bt & sqrtIV >=sqrtIVt & p.value<=alpha
    A = A[sel]
    B = B[sel]
    sqrtIV = sqrtIV[sel]
    p.value = p.value[sel]
    nc = 4
  }
  sel[is.na(sel)]=FALSE
  CM = subset(C,sel)
  m = data.frame(matrix(0,nrow = sum(sel), ncol=nc))
  if(nc==4) names(m) = c("A","B","sqrtIV", "p.value")
  else if(nc==10) names(m) = c("A","LA","UA","B","LB","UB","sqrtIV", "LsqrtIV","UsqrtIV", "p.value")
  if(sum(sel)>0) {
    if(nc==4) {
      m[,1] = A
      m[,2] = B
      m[,3] = sqrtIV
      m[,4] = p.value
    } else if (nc==10){
      m[,1:3] = A
      m[,4:6] = B
      m[,7:9] = sqrtIV
      m[,10] = p.value
    }
    row.names(m) = rownames(CM)
  }
  m = m[order(m[["sqrtIV"]],decreasing=TRUE),]
  print(m,...)
}

#' @rdname indicators
#' @param object An object of class 'indicators'.
#' @param ... Additional arguments for functions \code{print},\code{summary} or \code{plot}.
summary.indicators <- function (object, ...) {
  cat(paste("Result of 'indicators()':\n\n"))
  cat(paste("  Number of plots in target site group: ", sum(object$group.vec),"\n",sep=""))
  cat(paste("  Number of candidate species: ", length(object$candidates),"\n",sep=""))
  cat(paste("  Number of species used in combinations: ", length(object$candidates),"\n",sep=""))
  cat(paste("  Number of indicators (single species or species combinations) kept: ", nrow(object$C),"\n",sep=""))
  cat(paste("  Group coverage: ", coverage(object)*100,"%\n",sep=""))
  cat("\n")
}

#' @rdname indicators
#' @param type Statistic to plot. Accepted values are "IV" (indicator value), "sqrtIV" (square root of indicator value), "A", "LA", "UA", (positive predictive value and confidence limits), "B",  "LB", "UB" (sensitivity and confidence limits).
#' @param maxline Flag to indicate whether a line has to be drawn joining the maximum values for each order of combinations.
plot.indicators<-function(x, type="sqrtIV", maxline=TRUE, ...) {
  A = x$A
  B = x$B
  sqrtIV=x$sqrtIV
  order = rowSums(x$C)
  if(is.data.frame(A)) {
    if(type=="IV") val = sqrtIV[,1]^2
    else if(type=="sqrtIV") val = sqrtIV[,1]
    else if(type=="A") val = A[,1]	
    else if(type=="B") val = B[,1]	
    else if(type=="LA") val = A[,2]	
    else if(type=="UA") val = A[,3]	
    else if(type=="LB") val = B[,2]	
    else if(type=="UB") val = B[,3]	
    else if(type=="LsqrtIV") val = sqrtIV[,2]	
    else if(type=="UsqrtIV") val = sqrtIV[,3]	
  } else {
    if(type=="IV") val = sqrtIV^2
    else if(type=="sqrtIV") val = sqrtIV
    else if(type=="A") val = A	
    else if(type=="B") val = B	   	
  }
  plot(order,val, type="n", axes=FALSE, xlab="Order", ylab=type,...)	
  points(order,val, pch=1, cex=0.5)	
  axis(1, at = order, labels=order)
  axis(2)
  if(maxline) lines(1:max(order),tapply(val,order,max), col="gray")
}
