#' Strength of species site-group associations
#'
#' This function computes the strength of the association between a species vector and a vector of memberships to the site-group. Several statistics are possible, following De \enc{Cáceres}{Caceres} and Legendre (2009). The choice of the most appropriate statistic depends on the application. If required, bootstrap confidence interval bounds are also computed.
#' 
#' @param X Community data table (rows are sites, columns are species).
#' @param cluster A vector of numeric group labels for sites.
#' @param func The association index to be used. Accepted function names: "r", "r.g", "IndVal", "IndVal.g", "A", "A.g", "B", "cos", "cos.g", "r.ind", "r.ind.g", "s.ind", "s.ind.g" (lowercase values are also accepted).
#' @param group Group for which association values are to be computed. If group=NULL, association values will be computed for all groups.
#' @param nboot.ci Number of bootstrap samples. If \code{NULL} then confidence intervals are not estimated.
#' @param alpha.ci Error in confidence intervals.
#' @param c Vector of total abundance per site (used in functions "r.ind", "r.ind.g", "s.ind", "s.ind.g" only).
#' 
#' @details
#' This R function is applicable to both presence-absence and quantitative species data, depending on the values in the input matrix X.  Indices "r","r.g","cos" and "cos.g" correspond to the general correlation generalization ("r" is the point-biserial correlation coefficient), whereas "r.ind", "r.ind.g", "s.ind" and "s.ind.g" correspond to the individual-based generalization. Both approaches give the same results for binary (presence-absence) species data. Indicator value indices "IndVal" and "IndVal.g" are partial generalizations of the indices of the presence-absence indices. The "IndVal" index of \enc{Dufrêne}{Dufrene} & Legendre (1997) is called "IndVal.g" here and \code{strassoc} returns the square root of the original index. Indices "A","A.g" and "B" are the asymmetric quantities into which indicator values "IndVal" and "IndVal.g" can be decomposed. See De Cáceres and Legendre for details on the relationships between these indices and their usage. Bootstrap confidence intervals are computed using the simple percentile method (Manly 1997).
#' 
#' @return
#' Returns a matrix of association values, where species are in rows and groups are in columns. If \code{nboot.ci} is not null, then a list is returned, where 'stat' contains the matrix of association values, and the lower and upper confidence limits are given in supplementary matrices named 'lowerCI' and 'upperCI' respectively.
#' 
#' @note
#' This function gives the same association values as function \code{indval} in package "labdsv" when used setting \code{func="IndVal.g"}, excepting the fact that the square root IndVal values is returned instead of the original IndVal.
#' 
#' @references 
#' De \enc{Cáceres}{Caceres}, M. and Legendre, P. 2009. Associations between species and groups of sites: indices and statistical inference. Ecology 90(12): 3566-3574.
#' 
#' \enc{Dufrêne}{Dufrene}, M. and P. Legendre. 1997. Species assemblages and indicator species: The need for a flexible asymetrical approach. Ecological Monographs 67:345-366.
#' 
#' Manly, B. F. J. 1997. Randomization, bootstrap and Monte Carlo methods in biology. Chapman and Hall Texts in Statistical Science Series.
#' 
#' @author Miquel De \enc{Cáceres}{Caceres} Ainsa,  EMF-CREAF
#' 
#' @seealso \code{\link{signassoc}}, \code{\link{multipatt}}
#' 
#' @export
#'
#' @examples
#' ## Load species data
#' data(wetland) 
#' 
#' ## Create three clusters using
#' wetkm <- kmeans(wetland, centers=3)
#' 
#' ## Compute Dufrene and Legendre's IndVal
#' strassoc(wetland, wetkm$cluster, func="IndVal.g") 
#' 
#' ## Compute correlation, with bootstrap 95 percent confidence intervals
#' strassoc(wetland, wetkm$cluster, func="r.g", 
#'          nboot.ci =100) 
#'          
strassoc <- function(X, cluster,func="r",group=NULL, nboot.ci=NULL, alpha.ci=0.05, c=1) {

  # Check arguments
  func <- match.arg(func, c("r", "r.g", "indval", "IndVal", "indval.g", "IndVal.g","a", "A", "a.g", "A.g", "B", "b","cos", "cos.g", "r.ind", "r.ind.g", "s.ind", "s.ind.g"))

  #Turn into a matrix (if not)
  X = as.matrix(X)
  
IndVal1<-function(sav,gmv,group=NULL) {
   gmv= as.factor(gmv)
   means = vector("numeric",length(levels(gmv)))
   for(i in 1:length(levels(gmv))) {
   		means[i]=mean(sav[gmv==levels(gmv)[i]])
   		if(is.na(means[i])) means[i]=0
   }   
   if(is.null(group)) { #COMPUTE INDVAL1 FOR ALL GROUPS 
	   indvals = data.frame(matrix(0, nrow=length(levels(gmv)),ncol=3))
   	   row.names(indvals)=levels(gmv)
   	   names(indvals)=c("A.g","B","IndVal.g")
   	   for(i in 1:length(levels(gmv))) {
	  		indvals[i,1]= ifelse(sum(means)>0,means[i]/sum(means),0)
	  		indvals[i,2] = ifelse(sum(gmv==levels(gmv)[i])>0, sum((gmv==levels(gmv)[i]) & (sav>0))/sum(gmv==levels(gmv)[i]),0)
    		indvals[i,3]=sqrt(indvals[i,1]*indvals[i,2])
   	   }	 
   } else { #COMPUTE INDVAL1 FOR ONE GROUP
	   indvals = data.frame(matrix(0, nrow=1,ncol=3))
   	   row.names(indvals)=c(group)
   	   names(indvals)=c("A.g","B","IndVal.g")
   	   m = mean(sav[gmv==group])
   	   if(is.na(m)) m=0
   	   indvals[1,1]= ifelse(sum(means)>0,m/sum(means),0)
	   indvals[1,2] = ifelse(sum(gmv==group)>0, sum((gmv==group) & (sav>0))/sum(gmv==group),0)
       indvals[1,3]=sqrt(indvals[1,1]*indvals[1,2])
   }
   return (indvals)
}

IndVal2<-function(sav,gmv,group=NULL) {
   gmv= as.factor(gmv)
   sums = vector("numeric",length(levels(gmv))) 
   for(i in 1:length(levels(gmv))) {
	   sums[i]=sum(sav[gmv==levels(gmv)[i]])
   }
   
   if(is.null(group)) {   #COMPUTE INDVAL2 FOR ALL GROUPS 
   	   indvals = data.frame(matrix(0, nrow=length(levels(gmv)),ncol=3))
       row.names(indvals)=levels(gmv)
	   
	   for(i in 1:length(levels(gmv))) {
		  indvals[i,1]= ifelse(sum(sums)>0,sums[i]/sum(sums),0)
		  indvals[i,2] = ifelse(sum(gmv==levels(gmv)[i])>0, sum((gmv==levels(gmv)[i]) & (sav>0))/sum(gmv==levels(gmv)[i]),0)
		  indvals[i,3]=sqrt(indvals[i,1]*indvals[i,2])
	   }
   } else {  #COMPUTE INDVAL2 FOR ONE GROUP
   	   indvals = data.frame(matrix(0, nrow=1,ncol=3))
       row.names(indvals)=c(group)
	   sg=sum(sav[gmv==group])
	   indvals[1,1]= ifelse(sum(sums)>0,sg/sum(sums),0)
	   indvals[1,2] = ifelse(sum(gmv==group)>0, sum((gmv==group) & (sav>0))/sum(gmv==group),0)
	   indvals[1,3]=sqrt(indvals[1]*indvals[2])       
   }
   names(indvals)=c("A","B","IndVal")
   
   return (indvals)
}

A<-function(sav,gmv,group=NULL) {
   gmv= as.factor(gmv)
   sums = vector("numeric",length(levels(gmv))) 
   for(i in 1:length(levels(gmv))) {
	   sums[i]=sum(sav[gmv==levels(gmv)[i]])
   }
   if(is.null(group)) {   #COMPUTE A FOR ALL GROUPS 
   	   A = rep(0,length(sums))
	   for(i in 1:length(levels(gmv))) {
		  A[i]= ifelse(sum(sums)>0,sums[i]/sum(sums),0)
	   }
   } else {  #COMPUTE A FOR ONE GROUP
	   A= ifelse(sum(sums)>0,sum(sav[gmv==group])/sum(sums),0)
   }   
   return (A)
}
r.g<-function(sav,gmv, group=NULL) {   
   gmv= as.factor(gmv)
   npm = vector("numeric",length(levels(gmv)))
   k=length(levels(gmv))
   N=length(sav)
   nm = 0
   l2m = 0
   npmg=0
   C=1/k
   for(i in 1:k) {
	   Np=sum(gmv==levels(gmv)[i])
	   np=sum(sav*(gmv==levels(gmv)[i]))
	   n2p=sum((sav^2)*(gmv==levels(gmv)[i]))
	   npm[i] = (np/Np)*(N*C)
	   if(Np==0) npm[i]=0
	   nm=nm+npm[i]
	   a= (n2p/Np)*(N*C)
	   if(Np==0) a=0
	   l2m=l2m+a	   
	   if(!is.null(group)) if(levels(gmv)[i]==group) npmg=npm[i]
   }   
   Nm = k*N*C   
   s2=(Nm*l2m)-(nm^2)
   g2=C*(1-C)
   if(is.null(group)) {
	   r = vector("numeric",length(levels(gmv)))
	   for(i in 1:k) {
			num=npm[i]-(nm*C)
			r[i]=num/sqrt(s2*g2)
			if(num==0) r[i]=0
			if(is.nan(r[i])) cat(l2m, " " , nm, "\n")
	   }	 
   } else {
	   #print((npmg/k)+nm)   	 
	   num=npmg-(nm*C)
	   r=num/sqrt(s2*g2)	
	   if(num==0) r=0	   	   
   }
   return(r)
}

r.s<-function(sav,gmv,group=NULL) {   
   gmv= as.factor(gmv)
   if(is.null(group)) {
	   r = vector("numeric",length(levels(gmv)))   
	   for(i in 1:length(levels(gmv))) {
			r[i]=cor(sav,gmv==levels(gmv)[i])
			if(is.na(r[i])) r[i]=0
	   }	 
   } else {
   		r=cor(sav,ifelse(gmv==group,1,0))
   		if(is.na(r)) r=0 #Si la desv estandar es zero no es pot calcular i posem r=0
   }
   return(r)
}

r.ind.s<-function(sav,gmv, group=NULL, c=1) {   
   gmv= as.factor(gmv)
   N=length(sav)
   asp=sum(sav)
   if(is.null(group)) {
	   r = vector("numeric",length(levels(gmv)))
	   for(i in 1:length(levels(gmv))) {
	       ni=sum(gmv==levels(gmv)[i])
	       aspi=sum(sav*(gmv==levels(gmv)[i]))
			num=(N*aspi)-(asp*ni)
			s2 = N*c*asp - asp^2
			g2 = N*ni - ni^2
			r[i]=num/sqrt(s2*g2)
			if(num==0) r[i]=0
	   }	 
   } else {
   	   ni=sum(gmv==group)
	   aspi=sum(sav*(gmv==group))
	   num=(N*aspi)-(asp*ni)
	   s2 = N*c*asp - asp^2
	   g2 = N*ni - ni^2
	   r=num/sqrt(s2*g2)
	   if(num==0) r=0	   	   
   }
   return(r)
}

r.ind.g<-function(sav,gmv, group=NULL, c=1) {   
   gmv= as.factor(gmv)
   aspig = vector("numeric",length(levels(gmv)))
   k= length(levels(gmv))
   N=length(sav)
   aspg = 0
   nig=N/k
   for(i in 1:k) {
	   aspig[i]=(N/k)*(sum(sav*(gmv==levels(gmv)[i]))/sum(gmv==levels(gmv)[i]))
	   aspg = aspg + aspig[i]
	   if(!is.null(group)) {
	   	   if(group==levels(gmv)[i]) igroup = i
	   	}
   }
   
   if(is.null(group)) {
	   r = vector("numeric",length(levels(gmv)))
	   for(i in 1:k) {
			num=(N*aspig[i])-(aspg*nig)
			s2 = N*c*aspg - aspg^2
			g2 = N*nig - nig^2
			r[i]=num/sqrt(s2*g2)
			if(num==0) r[i]=0
	   }	 
   } else {
		num=(N*aspig[igroup])-(aspg*nig)
	   s2 = N*c*aspg - aspg^2
	   g2 = N*nig - nig^2
	   r=num/sqrt(s2*g2)
	   if(num==0) r=0	   	   
   }
   return(r)
}

s.ind.s<-function(sav,gmv, group=NULL, c=1) {   
   gmv= as.factor(gmv)
   N=length(sav)
   asp=sum(sav)
   if(is.null(group)) {
	   r = vector("numeric",length(levels(gmv)))
	   for(i in 1:length(levels(gmv))) {
	       ni=sum(gmv==levels(gmv)[i])
	       aspi=sum(sav*(gmv==levels(gmv)[i]))
			r[i]=aspi/sqrt(asp*c*ni)
			if(aspi==0) r[i]=0
	   }	 
   } else {
   	   ni=sum(gmv==group)
	   aspi=sum(sav*(gmv==group))
		r=aspi/sqrt(asp*c*ni)
	   if(aspi==0) r=0	   	   
   }
   return(r)
}

s.ind.g<-function(sav,gmv, group=NULL, c=1) {   
   gmv= as.factor(gmv)
   N=length(sav)
   k=length(levels(gmv))
   asp=0
   for(i in 1:k) {
  	  asp = asp + sum(sav*(gmv==levels(gmv)[i]))/sum(gmv==levels(gmv)[i])
   }	 
   if(is.null(group)) {
	   r = vector("numeric",length(levels(gmv)))
	   for(i in 1:k) {
	       ni=sum(gmv==levels(gmv)[i])
	       aspi=sum(sav*(gmv==levels(gmv)[i]))
			r[i]=sqrt(((aspi/ni)*aspi)/(asp*c*ni))
			if(aspi==0) r[i]=0
	   }	 
   } else {
   	   ni=sum(gmv==group)
	   aspi=sum(sav*(gmv==group))
#	   r=sqrt(((aspi/ni)*aspi)/(asp*c*ni))
	   r=sqrt(((aspi/ni))/(asp*c))
	   if(aspi==0) r=0	   	   
   }
   return(r)
}

cos.g<-function(sav,gmv, group=NULL) {   
   gmv= as.factor(gmv)
   s = 0
   k=length(levels(gmv))
   if(is.null(group)) {
	   n = vector("numeric",length(levels(gmv)))
	   a = vector("numeric",length(levels(gmv)))
  	   r = vector("numeric",length(levels(gmv)))
   	   for(i in 1:k) {
      		n[i]=sum(gmv==levels(gmv)[i])
	   		a[i]=sum(sav*(gmv==levels(gmv)[i]))
	   		a2=sum((sav^2)*(gmv==levels(gmv)[i]))
	   		s = s+(a2/n[i])
   		}
	   for(i in 1:k) r[i]=(a[i]/n[i])/sqrt(s)
   	} else {
   	   for(i in 1:k) {
      		n =sum(gmv==levels(gmv)[i])
	   		a2=sum((sav^2)*(gmv==levels(gmv)[i]))
	   		s = s+(a2/n)
   		}
	   r=(sum(sav*(gmv==group))/sum(gmv==group))/sqrt(s)
   	}
   return(r)
}


cos.s<-function(sav, gmv, group=NULL) {
   gmv= as.factor(gmv)
   r = vector("numeric",length(levels(gmv)))
   x1 <- as.numeric(sav)
   l1=sqrt(sum(x1*x1))
   if(is.null(group)) {
   		k=length(levels(gmv))
   		for(i in 1:k) {
      		x2 <- ifelse(gmv==levels(gmv)[i],1,0)
	  		l2=sqrt(sum(x2*x2))
	  		r[i] = (sum(x1*x2)/(l1*l2))    
	  		if(is.na(r[i])) r[i]=0 
   		}
   	} else {
   		x2 <- ifelse(gmv==group,1,0)
	  	l2=sqrt(sum(x2*x2))
	  	r = (sum(x1*x2)/(l1*l2))    
	  	if(is.na(r)) r=0 
   	}
   return(r)
}

  ass.func<-function(sav,gmv, func="r", group=NULL) {
	if(func=="r") a = r.s(sav,gmv, group)
	else if(func=="r.g") a= r.g(sav,gmv,group)
	else if(func=="cos") a = cos.s(sav,gmv,group)
	else if(func=="cos.g") a = cos.g(sav,gmv,group)
	else if(func=="r.ind") a= r.ind.s(sav,gmv,group,c)
	else if(func=="r.ind.g") a= r.ind.g(sav,gmv,group,c)
	else if(func=="s.ind") a= s.ind.s(sav,gmv,group,c)
	else if(func=="s.ind.g") a= s.ind.g(sav,gmv,group,c)
	else if(func=="IndVal.g" || func=="indval.g") a = IndVal1(sav,gmv,group)[,3]
	else if(func=="IndVal" || func == "indval") a = IndVal2(sav,gmv,group)[,3]
	else if(func=="A.g" || func == "a.g") a = IndVal1(sav,gmv,group)[,1]
	else if(func=="A" || func == "a") a = A(sav,gmv,group)
	else if(func=="B" || func == "b") a = IndVal2(sav,gmv,group)[,2]
	return (a)
  }	
  
  if(sum(is.na(cluster))>0) stop("Cannot deal with NA values. Remove and run again.")
  if(sum(is.na(X))>0) stop("Cannot deal with NA values. Remove and run again.")
  
  cluster= as.factor(cluster)
  if(is.matrix(X)| is.data.frame(X)) {
  	 nsps = dim(X)[2]
	 nsites = dim(X)[1]
  }
  else {
  	 nsps = 1
  	 nsites = length(X)
  }
  ngroups = length(levels(cluster))	
  if(!is.null(group)) ngroups=1
  dm = matrix(0,nsps,ngroups)
  
  
  rownames(dm)=colnames(X)
  if(is.null(group)) colnames(dm)=levels(cluster)
  else colnames(dm)=group
  
	# Compute diagnostic values
	for(i in 1:nsps) {	
		dm[i,] = ass.func(X[,i], cluster,func,group)
	}
	
	if(!is.null(nboot.ci)) {
		dmb = array(0,dim=c(nboot.ci,nsps,ngroups))
		dmlower = matrix(0,nsps,ngroups)
		dmupper = matrix(0,nsps,ngroups)
		for(b in 1:nboot.ci) {
		  bi = sample(nsites,replace=TRUE)
			clusterb = cluster[bi]
			for(i in 1:nsps) {	
				savb = X[bi,i]
				dmb[b,i,] = ass.func(savb, clusterb,func,group)
			}
		}

		for(i in 1:nsps) {	
			for(k in 1:ngroups) {	
				   sdmb = sort(dmb[,i,k])
				   dmlower[i,k]=sdmb[(alpha.ci/2.0)*nboot.ci]
				   dmupper[i,k]=sdmb[(1-(alpha.ci/2.0))*nboot.ci]
			}
		}
    rownames(dmlower)=colnames(X)
	   if(is.null(group)) colnames(dmlower)=levels(cluster)
	   else colnames(dmlower)=group
  	 rownames(dmupper)=colnames(X)
	   if(is.null(group)) colnames(dmupper)=levels(cluster)
	   else colnames(dmupper)=group
	   
	   return(list(stat=dm,lowerCI=dmlower,upperCI=dmupper))
    }
   return(dm)
 }

