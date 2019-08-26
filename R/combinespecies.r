combinespecies<-function(X, min.order = 1, max.order = 3, min.occ = 0, FUN = min, verbose=FALSE, add.names=TRUE, ...) {

  if(sum(is.na(X))>0) stop("Cannot deal with NA values. Remove and run again.")
  if(min.order<1) stop("min.order must be at least 1.")
  if(min.order>max.order) stop("min.order must be equal or larger than max.order.")
  nsites = nrow(X)
  
  #Get species names
  spplist = names(X)
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
  
  C = data.frame(matrix(0,nrow=nc,ncol=nspecies))
  names(C)<-spplist
  XC = data.frame(matrix(0, nrow=nrow(X), ncol=nc))
  row.names(XC) = row.names(X)
  for(i in 1:nc) {
    spvec = as.numeric(comblistDef[[i]])
    XC[,i] <-apply(X[,spvec, drop=FALSE],1,FUN=FUN,...)
    C[i,spvec] = 1
  }
  if(add.names) {
    for(i in 1:ncol(XC)) names(XC)[i] = paste(spplist[C[i,]==1], collapse="+")
  } else names(XC) = 1:ncol(XC)
  row.names(C) = names(XC)
  return(list(XC = XC, C=C))
}