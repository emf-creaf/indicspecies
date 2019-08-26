combinespecies<-function(X, min.order = 1, max.order = 3, min.occ = 1, FUN = min, verbose=FALSE, add.names=TRUE, ...) {

  #Turn into a matrix (if not)
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