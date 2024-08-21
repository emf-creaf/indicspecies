#' Multi-level pattern analysis
#' 
#' This function studies the association between species patterns and combinations of groups of sites.
#'
#' @param x Community data table
#' @param cluster A vector representing a partition of sites
#' @param func Species-site group association function. Four values are accepted \code{"IndVal"}, \code{"IndVal.g"}, \code{"r"} and \code{"r.g"} (lowercase values are also accepted).
#' @param duleg If TRUE, site group combinations are not considered, only the original site groups, like in Dufrêne & Legendre (1997). Internally, \code{duleg = TRUE} equals \code{max.order = 1}.
#' @param restcomb A vector of integer values used to restrict the combinations of site groups to those with ecological sense according to the analyst. The default \code{NULL} indicates that all combinations are used. If \code{duleg=TRUE} this argument is ignored.
#' @param min.order An integer indicating the minimum order of site group combinations (by default \code{max.order=1} for singletons). Cannot be larger than \code{max.order}.
#' @param max.order An integer indicating the maximum order of site group combinations to be considered: \code{max.order=1} for singletons, \code{max.order=2} for pairs, \code{max.order=3} for triplets... As \code{restcomb}, this parameter provide a way to restrict the site group combinations that make ecological sense. By default all possible site group combinations are considered. If \code{max.order=1} then the function will behave as if \code{duleg=TRUE}.
#' @param control A list of control values describing properties of the permutation design, as returned by a call to \code{\link[permute]{how}}.
#' @param permutations A custom matrix of permutations, to be used if \code{control = NULL}, with permutations in rows and site indices in columns
#' @param print.perm If TRUE, prints permutation numbers after each set of 100 permutations.
#' 
#' @details
#' This function creates combinations of the input clusters and compares each combination with the species in the input matrix x. For each species it chooses the combination with a highest association value. Best matching patterns are tested for statistical significance of the associations. Four association indices are possible (some less than for \code{\link{strassoc}}): "IndVal", "IndVal.g", "r" and "r.g". Indicator value indices will return the pattern that better matches the species observed pattern, whereas correlation indices will return the pattern that creates a highest inside/outside difference. Details are given in De \enc{Cáceres}{Caceres} et al. (2010). The user can restrict the combinations in three ways: (1) by using \code{duleg=TRUE}, which leads to consider single site-groups only; (2) by setting the minimum and maximum order of combinations using \code{min.order} and \code{max.order}; or (3) by using \code{restcomb} to restrict combinations at will. In order to carry out the third way, values in \code{restcomb} must be the indices of combinations that appear in the column \code{index} of object \code{sign} (see below). 
#' 
#' Complex permutation designs are allowed through the function \code{\link[permute]{how}} from package "permute". If those are not enough, the user can set \code{control = NULL} and specify a custom matrix of permutations to test with parameter \code{permutations}.
#' 
#'
#' @return
#' An object of class \code{multipatt} with:
#' \item{func}{The name of the function used.}
#' \item{comb}{A matrix describing all the combinations studied.}
#' \item{str}{A matrix the association strength for all combinations studied.}
#' \item{A}{If \code{func = "IndVal"} (or \code{func = "IndVal.g"}) a matrix whose values are the "A" (or "A.g") component of indicator values. Otherwise this element is left as \code{NULL}.}
#' \item{B}{If \code{func = "IndVal"} (or \code{func = "IndVal.g"}) a matrix whose values are the "B" component of indicator values. Otherwise this element is left as \code{NULL}.}
#' \item{sign}{Data table with results of the best matching pattern, the association value and the degree of statistical significance of the association (i.e. p-values from permutation test). Note that p-values are not corrected for multiple testing.}
#' 
#' @note
#' This function gives the same results as function \code{indval} in package "labdsv" when used setting \code{func="IndVal.g"} and \code{duleg=TRUE}, excepting the fact that the square root IndVal values is returned instead of the original IndVal.
#' 
#' @references 
#' De \enc{Cáceres}{Caceres}, M. and Legendre, P. 2009. Associations between species and groups of sites: indices and statistical inference. Ecology 90(12): 3566-3574.
#' 
#' De \enc{Cáceres}{Caceres}, M., Legendre, P., Moretti, M. 2010. Improving indicator species analysis by combining groups of sites. Oikos 119(10): 1674-1684.
#' 
#' \enc{Dufrêne}{Dufrene}, M. and P. Legendre. 1997. Species assemblages and indicator species: The need for a flexible asymetrical approach. Ecological Monographs 67:345-366.
#' 
#' @author 
#' Miquel De \enc{Cáceres}{Caceres} Ainsa, EMF-CREAF
#' 
#' Florian Jansen, Institute of Botany and Landscape Ecology, Ernst-Moritz-Arndt-University
#' 
#' @seealso \code{\link{summary.multipatt}}, \code{\link{strassoc}}, \code{\link{signassoc}}, \code{\link[permute]{how}}
#' 
#' @export
#'
#' @examples
#' library(stats)
#' 
#' data(wetland) ## Loads species data
#' 
#' wetkm <- kmeans(wetland, centers=3) ## Creates three clusters using kmeans
#' 
#' ## Runs the combination analysis using IndVal.g as statistic
#' wetpt <- multipatt(wetland, wetkm$cluster, control = how(nperm=999)) 
#' 
#' ## Lists those species with significant association to one combination
#' summary(wetpt) 
#' 
#' ## Lists those species with significant association to one combination, 
#' ## including indval components.
#' summary(wetpt, indvalcomp=TRUE) 
#' 
multipatt <- function (x, cluster, func = "IndVal.g", duleg = FALSE, restcomb=NULL, 
                       min.order=1, max.order=NULL, 
                       control = how(), permutations = NULL, print.perm = FALSE)                                                                                              
{
	                                                                                                                              

# Matrix of possible cluster combinations
cl.comb <- function(clnames, min.order, max.order) {
	  k <- length(clnames)
    ep <- NULL
    names.ep <- NULL
    for(j in max(min.order,1):min(max.order, k)) {
      nco <- choose(k,j)
      co <- combn(k,j)
      epn <- matrix(0,ncol=nco,nrow=k)
      for(i in 1:ncol(co)) {
        epn[co[,i],i] <- 1
        if(is.null(names.ep)) names.ep <- paste(clnames[co[,i]], collapse = "+")
        else names.ep <- c(names.ep, paste(clnames[co[,i]], collapse = "+"))
      }
      if(is.null(ep)) ep <- epn
      else ep <- cbind(ep,epn)
    }
    colnames(ep) <- names.ep
    return(ep)
}

# Correlation measures for combinations
rcomb <- function(x, memb, comb, min.order, max.order, mode="group", restcomb=NULL) {
  k = ncol(memb)
  nsps = ncol(x) #Number of species
  N = dim(comb)[1]	#Number of sites
  # ni = diag(t(memb) %*% memb)
  ni = colSums(memb*memb)
  tx <- t(x)
  aisp = (tx %*% memb)
  lisp = (tx^2 %*% memb)
  
  if(mode=="site") {
    lspK = rowSums(lisp)
    aspK = rowSums(aisp)		
    aspC = (tx %*% comb)
    # nC = diag(t(comb) %*% comb) too much memory demanded!
    nC = colSums(comb*comb)
  } else if(mode=="group") {
    aispni=sweep(aisp,2,ni,"/")
    lispni=sweep(lisp,2,ni,"/")
    #Corrected sum and length of species vectors
    lspK = (N/k) * rowSums(lispni)
    aspK = (N/k) * rowSums(aispni)
    if(max.order==1) {
      aspC = (N/k)*aispni
      nC = rep(N/k,k)
    } else {
      #Corrected sum of species values in combinations
      aspC = matrix(0,nrow=nsps,ncol=ncol(comb))
      #Corrected size of cluster combinations
      nC = vector(mode="numeric",length=ncol(comb))
      
      cnt = 1
      for(level in min.order:max.order) {
        co = combn(1:k,level)
        for(j in 1:ncol(co)) {
          if(nrow(co)>1) aspC[,cnt] = rowSums(aispni[,co[,j]])
          else aspC[,cnt] = aispni[,co[,j]]
          nC[cnt] = length(co[,j])
          cnt= cnt+1
        }
      }
      aspC = (N/k)*aspC
      nC = (N/k)*nC
    }
  }
  #Compute index
  num = N*aspC - aspK%o%nC
  den = sqrt(((N*lspK)-aspK^2)%o%(N*nC-(nC^2)))
  str=num/den
  colnames(str) <- colnames(comb)[1:ncol(str)]
  #Remove site group combinations that are not to be compared
  if(!is.null(restcomb)) {
    if(sum(restcomb %in% (1:ncol(str)))!=length(restcomb)) 
      stop(paste("One or more indices in 'restcomb' are out of range [1, ",ncol(str),"]",sep=""))      
    str <- str[,restcomb, drop = FALSE]
  }
  return(str)
}

# IndVal for combinations
indvalcomb <- function(x, memb, comb, min.order, max.order, mode = "group", restcomb=NULL, indvalcomp=FALSE) {
  k <- ncol(memb)
  tx <- t(x)
  aisp = tx %*% comb
  dx <- dim(tx)
  nisp <- matrix(as.logical(tx),nrow=dx[1],ncol=dx[2]) %*% comb
  # ni = diag(t(comb) %*% comb)
  ni = colSums(comb*comb)
  nispni = sweep(nisp, 2, ni, "/")   
  if (mode == "site") A = sweep(aisp, 1, colSums(x), "/")  
  else {
    aispni = sweep(tx %*% memb, 2, colSums(memb), "/")
    asp = rowSums(aispni)
    if(max.order==1) {
      A <- sweep(aispni, 1, asp, "/")
    } else {
      s <- NULL 
      for(j in min.order:max.order) {
        if(j == 1) {
          s <- aispni
        } else {
          co <- combn(k,j)
          sn <- apply(co, 2, function(x) rowSums(aispni[,x]))
          if(is.null(s)) s <- sn
          else s <- cbind(s, sn)          
        }
      }
	    A = sweep(s, 1, asp, "/")
  	} 
  }
  iv = sqrt(A * nispni)
  colnames(iv) <- colnames(comb)
  colnames(A) <- colnames(comb)
  colnames(nispni) <- colnames(comb)
  rownames(A) <- rownames(iv)
  rownames(nispni) <- rownames(iv)

  #Remove site group combinations that are not to be compared
  if(!is.null(restcomb)) {
    if(sum(restcomb %in% (1:ncol(iv)))!=length(restcomb)) 
      stop(paste("One or more indices in 'restcomb' are out of range [1, ",ncol(iv),"]",sep=""))
    iv = iv[,restcomb, drop = FALSE]
    A = A[,restcomb, drop = FALSE]
    nispni = nispni[,restcomb, drop = FALSE]
  }
  if(!indvalcomp) return(iv)
  else return(list(A=A,B=nispni, iv=iv))
}


  #function multipatt starts here	
  x <- as.matrix(x)                                                                                                                                
  vegnames <- colnames(x)
  nsps = ncol(x)
  clnames = levels(as.factor(cluster))
  k = length(clnames)
  
  if(!is.null(control)){
    nperm = control$nperm
  } else if (!is.null(permutations)){
    nperm = nrow(permutations)
  } else {
    stop("You must control permutations with how() or supply a matrix of permutations")
  }
  
  # Check parameters
  func= match.arg(func, c("r","r.g","IndVal.g","IndVal", "indval", "indval.g"))
  if(k<2) stop("At least two clusters are required.")
  if(sum(is.na(cluster))>0) stop("Cannot deal with NA values. Remove and run again.")
  if(sum(is.na(x))>0) stop("Cannot deal with NA values. Remove and run again.")
  if(!is.null(min.order)) {
    if(!is.numeric(min.order)) stop("Parameter min.order must be an integer.")
    if(func=="IndVal.g" || func=="indval.g" || func=="IndVal"|| func=="indval") min.order = min(max(round(min.order),1),k)
    else min.order = min(max(round(min.order),1),k-1)
  } else {min.order=1}
  if(is.null(max.order)) {
    if(func=="IndVal.g" || func=="indval.g" || func=="IndVal"|| func=="indval") max.order=k
    else max.order = k-1
  } else {
    if(!is.numeric(max.order)) stop("Parameter max.order must be an integer.")
    if(func=="IndVal.g" || func=="indval.g" || func=="IndVal"|| func=="indval") max.order = min(max(round(max.order),1),k)
    else max.order = min(max(round(max.order),1),k-1)
    if(max.order < min.order) stop("The value of 'max.order' cannot be smaller than that of 'min.order'.")
  }
  if(!is.null(restcomb)) { 
    restcomb = as.integer(restcomb)
  }

  # creates combinations from clusters
  if(duleg) max.order = 1 #duleg = TRUE overrides max.order

  # possible combinations (can also be used for permutations)
  combin <- cl.comb(clnames, min.order, max.order)	

  # discard combinations to discard that cannot be calculated because of order limitation
  restcomb = restcomb[restcomb<=ncol(combin)] 

  #Builds the plot membership matrix corresponding to combinations
  clind = apply(sapply(clnames,"==",cluster),1,which)
  comb <- combin[clind,]

  #Original Membership matrix 
  memb<-diag(1,k,k)[clind,]

  # Computes association strength for each group
  A = NULL
  B = NULL
  if(func=="r") str = rcomb(x, memb, comb, min.order, max.order, mode = "site", restcomb = restcomb)
  else if(func=="r.g") str = rcomb(x, memb,comb, min.order, max.order, mode = "group", restcomb = restcomb)
  else if(func=="IndVal"|| func=="indval") {
  	  IndVal = 	indvalcomb(x, memb, comb, min.order, max.order, mode = "site", restcomb = restcomb, indvalcomp=TRUE)
  	  str = IndVal$iv
  	  A = IndVal$A
  	  B = IndVal$B
  	}
  else if(func=="IndVal.g"||func=="indval.g") {
  	  IndVal = 	indvalcomb(x,  memb, comb, min.order, max.order, mode = "group", restcomb = restcomb, indvalcomp=TRUE)
  	  str = IndVal$iv
  	  A = IndVal$A
  	  B = IndVal$B
  	}

  # Maximum association strength
  maxstr = apply(str,1,max) 
  wmax <- max.col(str)
  #prepares matrix of results
  if(!is.null(restcomb))  m <- as.data.frame(t(combin[,restcomb, drop = FALSE][,wmax]))
  else  m <- as.data.frame(t(combin[,wmax]))
  row.names(m) <- vegnames
  names(m) <- sapply(clnames, function(x) paste("s", x, sep='.'))
  m$index <- wmax
  m$stat <- apply(str,1,max)

  #Perform permutations and compute p-values
  pv <- 1
  for (p in 1:nperm) {
    if(!is.null(control)){
      pInd = shuffle(length(cluster), control=control)
    } else {
      pInd = permutations[p,]
    }
    tmpclind = clind[pInd]
    combp = combin[tmpclind,]
    membp = diag(1,k,k)[tmpclind,]
    tmpstr <- switch(func,
                     r   = rcomb(x, membp, combp, min.order, max.order, mode = "site", restcomb = restcomb),
                     r.g = rcomb(x, membp, combp, min.order, max.order, mode = "group", restcomb = restcomb),
                     indval = indvalcomb(x, membp, combp, min.order, max.order, mode = "site", restcomb = restcomb),
                     IndVal = indvalcomb(x, membp, combp, min.order, max.order, mode = "site", restcomb = restcomb),
                     indval.g= indvalcomb(x, membp, combp, min.order, max.order, mode = "group", restcomb = restcomb),
                     IndVal.g= indvalcomb(x, membp, combp, min.order, max.order, mode = "group", restcomb = restcomb)
    )
    tmpmaxstr <- vector(length=nrow(tmpstr))
    for(i in 1:nrow(tmpstr)) tmpmaxstr[i] <- max(tmpstr[i,])	# apply is more slowly in this case
    pv = pv + (tmpmaxstr >= m$stat)
  }
  
  m$p.value <- pv/(1 + nperm)
  #Put NA for the p-value of species whose maximum associated combination is the set of all combinations
  m$p.value[m$index == (2^k-1)] <- NA

  if(!is.null(restcomb))  comb<-comb[,restcomb, drop = FALSE]
  
  a = list(call=match.call(), func = func, cluster = cluster, comb = comb, str = str, A=A, B=B, sign = m)
  class(a) = "multipatt"
  return(a)
}


