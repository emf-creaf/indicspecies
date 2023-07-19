#' Predicts site group from indicators
#' 
#' Function \code{predict.indicators} takes an object of class \code{\link{indicators}} and determines the probability of the indicated site group given a community data set. If no new data set is provided, the function can calculate the probabilities corresponding to the original sites used to build the \code{\link{indicators}} object.
#'
#' @param object An object of class 'indicators'.
#' @param newdata A community data table (with sites in rows and species in columns) for which predictions are needed. This table can contain either presence-absence or abundance data, but only presence-absence information is used for the prediction. If \code{NULL}, then the original data set used to derive the \code{\link{indicators}} object is used as data.
#' @param cv A boolean flag to indicate that probabilities should be calculated using leave-one-out cross validation (i.e recalculating positive predictive value of indicators after excluding the target site).
#' @param ... In function \code{predict}, additional arguments not used (included for compatibility with \code{\link{predict}}). 
#'
#' @details
#' Function \code{\link{indicators}} explores the indicator value of the simultaneous occurrence of sets of species (i.e. species combinations). The method is described in De \enc{Cáceres}{Caceres} et al. (2012) and is a generalization of the Indicator Value method of \enc{Dufrêne}{Dufrene} & Legendre (1997). The current function \code{predict.indicators} is used to predict the indicated site group from the composition of a new set of observations. For communities where one or more of the indicator species combinations are found, the function returns the probability associated to the indicator that has the highest positive predictive value (if confidence intervals are available, the maximum value is calculated across the lower bounds of the confidence interval). For communities where none of the indicator species combinations is found, the function returns zeroes. 
#' If \code{newdata = NULL}, the function can be used to evaluate the predictive power of a set of indicators in a cross-validated fashion. For each site in the data set, recalculates the predictive value of indicators after excluding the information of the site, and then evaluates the probability of the site group.
#' 
#' @return
#' If confidence intervals are available in \code{x}, function \code{predict.indicators} returns a matrix where communities are in rows and there are three columns, correspoinding to the probability of the indicated site group along with the confidence interval. If confidence intervals are not available in \code{x}, or if \code{cv = TRUE}, then \code{predict.indicators} returns a single vector with the probability of the indicated site group for each community. 
#' 
#' @references 
#' De \enc{Cáceres}{Caceres}, M., Legendre, P., Wiser, S.K. and Brotons, L. 2012. Using species combinations in indicator analyses. Methods in Ecology and Evolution 3(6): 973-982.
#' 
#' \enc{Dufrêne}{Dufrene}, M. and P. Legendre. 1997. Species assemblages and indicator species: The need for a flexible asymetrical approach. Ecological Monographs 67:345-366.
#' 
#' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, EMF-CREAF
#' 
#' @seealso \code{\link{indicators}}, \code{\link{pruneindicators}} \code{\link{coverage}}, \code{\link{multipatt}}, \code{\link{strassoc}}, \code{\link{signassoc}}
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
#' 
#' ## Run indicator analysis with species combinations for the first group
#' sc <- indicators(X=wetland, cluster=wetkm$cluster, group=1, verbose=TRUE, At=0.5, Bt=0.2)
#' 
#' ## Use the indicators to make predictions of the probability of group #1
#' ## Normally an independent data set should be used, because 'wetland' was used to derive
#' ## indicators. The same would be obtained calling 'predict(sc)' without further arguments.
#' p <- predict(sc, wetland)
#' 
#' ## Calculate cross-validated probabilities (recalculates 'A' statistics once for each site 
#' ## after excluding it, and then calls predict.indicators for that site)
#' pcv <- predict(sc, cv = TRUE)
#' 
#' ## Show original membership to group 1 along with (resubstitution) predicted probabilities  
#' ## and cross-validated probabilities. Cross-validated probabilities can be lower for sites
#' ## originally belonging to the target site group and higher for other sites.
#' data.frame(Group1 = as.numeric(wetkm$cluster==1), Prob = p, Prob_CV = pcv)
#' 
predict.indicators<-function(object, newdata = NULL, cv = FALSE, ...) {
  C <- as.matrix(object$C)
  nc<-nrow(C)
  taxnames <-colnames(C)
  if(cv) {
    X <- object$X
    group.vec <- object$group.vec
    ng = sum(group.vec)
    func <- object$func
    cluster <- object$cluster
    nr = nrow(X)
    Astat<-matrix(0, nrow=nr, ncol=nc)
    for(r in 1:nr) {
      for(i in 1:nc) {
        spvec = colnames(C)[C[i,]==1]
        found<-sum(X[r,spvec, drop=FALSE]>0)==sum(C[i,]) #Is the species combination there? (for each site)
        if(found) {
          # Recalculate A excluding the row
          sc.ab <-apply(X[-r,spvec, drop=FALSE],1,min)
          scg = sc.ab[group.vec[-r]]
          if(func=="IndVal.g") {
            mg = (sum(scg)/ng)
            Astat[r,i] = mg/sum(tapply(sc.ab,cluster[-r], "mean"))
          } else {
            Astat[r,i] = sum(scg)/sum(sc.ab)
          }
        }
      }
    }
    p = apply(Astat,1,max)
    names(p)<-rownames(X)
  } else {
    if(is.null(newdata)) newdata = object$X
    if(sum(taxnames %in% colnames(newdata))<length(taxnames)) stop("Not all taxon names that form indicators could be found in compositional table.")
    postaxa <- numeric(length(taxnames))
    for(i in 1:length(taxnames)) postaxa[i] = which(colnames(newdata)==taxnames[i])
    stat<-matrix(0, nrow=nrow(newdata), ncol=nc)
    if(length(dim(object$A))==2) {
      lowerCI<-matrix(0, nrow=nrow(newdata), ncol=nc)
      upperCI<-matrix(0, nrow=nrow(newdata), ncol=nc)    
    }
    for(i in 1:nc) {
      ci <- as.numeric(C[i,])
      found<-rowSums(newdata[,postaxa[ci==1], drop=FALSE]>0)==sum(ci) #Is the species combination there? (for each site)
      if(length(dim(object$A))==2) {
        stat[,i]<- found*object$A$stat[i]
        lowerCI[,i]<- found*object$A$lowerCI[i]
        upperCI[,i]<- found*object$A$upperCI[i]
      }
      else {
        stat[,i]<- found*object$A[i]
      }
    }
    ## Returns the maximum probability among all combinations explored
    ## If confidence intervals are available
    if(length(dim(object$A))==2) {
      wm = apply(lowerCI,1,which.max)
      p = matrix(NA,nrow=nrow(newdata), ncol=3)
      for(i in 1:length(wm)) {
        p[i,] = c(stat[i,wm[i]],lowerCI[i,wm[i]], upperCI[i,wm[i]])
      }
      rownames(p)<-rownames(newdata)
      colnames(p)<-c("Prob.", "lowerCI", "upperCI")
    } else {
      p = apply(stat,1,max)
      names(p)<-rownames(newdata)
    }
  }
  return(p)
}