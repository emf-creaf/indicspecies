#' Summarizing multi-level pattern analysis
#' 
#' This is a \code{summary} function for \code{multipatt} objects.
#'
#' @param object An object returned by function \code{\link{multipatt}}.
#' @param alpha Significance level for selecting species in the summary.
#' @param minstat Minimum (absolute) value of the statistic for selecting species in the summary.
#' @param At Minimum value of positive predictive value (A) for selecting indicators. This argument is effective only if \code{multipatt} function was executed using \code{func = "IndVal"} or \code{func = "IndVal.g"}.
#' @param Bt Minimum value for sensitivity (B) for selecting indicators. This argument is effective only if \code{multipatt} function was executed using \code{func = "IndVal"} or \code{func = "IndVal.g"}.
#' @param indvalcomp Logical flag used to report components "A" (or "A.g") and "B" of indicator values. This argument is effective only if \code{multipatt} function was executed using \code{func = "IndVal"} or \code{func = "IndVal.g"}.
#' @param ... Additional arguments of the \code{summary} function.
#'
#' @details
#' The summary function for \code{multipatt} objects presents the results of the analysis by grouping species following their associated site group combination. Only significantly associated species are listed for each combination. The user can limit the species shown by setting a significance level and the minimum value for the strength of the association.
#' 
#' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, EMF-CREAF
#' 
#' @seealso \code{\link{strassoc}}, \code{\link{signassoc}}, \code{\link{multipatt}}
#' 
#' @export
#'
#' @examples
#' library(stats)
#' 
#' ## Loads species data
#' data(wetland) 
#' 
#' ## Creates three clusters using kmeans
#' wetkm <- kmeans(wetland, centers=3) 
#' 
#' ## Runs the combination analysis using IndVal.g as statistic
#' wetpt <- multipatt(wetland, wetkm$cluster) 
#' 
#' ## Lists those species with significant association to one combination
#' summary(wetpt) 
#' 
#' ## Lists those species with significant association to one combination, 
#' ## including indval components
#' summary(wetpt, indvalcomp=TRUE) 
#' 
summary.multipatt <- function (object, alpha = 0.05, minstat = NULL, At = NULL, Bt=NULL, 
                               indvalcomp=FALSE,...) {
  x <- object
  ncomb = ncol(x$str)
  ncolsign = ncol(x$sign)
  nsps = nrow(x$sign)
  cat("\n Multilevel pattern analysis")
  cat("\n ---------------------------\n")
  cat("\n Association function:", x$func)
  cat("\n Significance level (alpha):", alpha)
  if(!is.null(minstat)) cat("\n Minimum statistic value (minstat):", minstat)
  if(x$func=="IndVal" || x$func=="IndVal.g") {
    if(!is.null(At)) cat("\n Minimum positive predictive value (At):", At)
    if(!is.null(Bt)) cat("\n Minimum sensitivity (Bt):", Bt)
  }
  cat("\n\n Total number of species:", nsps)
  sel = !is.na(x$sign$p.value) & x$sign$p.value <= alpha
  if(!is.null(minstat)) sel = sel & (abs(x$sign$stat) >= minstat)
  if(!is.null(Bt) && !is.null(x$B)) {
    for(i in 1:nrow(x$sign)) sel[i] = sel[i] && (x$B[i,x$sign$index[i]]>=Bt)
  }
  if(!is.null(At) && !is.null(x$A)) {
    for(i in 1:nrow(x$sign)) sel[i] = sel[i] && (x$A[i,x$sign$index[i]]>=At)
  }
  a = x$sign[sel, ] #Selected matrix
  cat("\n Selected number of species:", nrow(a), "\n")
  cols = (ncolsign - 1):ncolsign
  
  #Collate indval components if required and possible
  if(indvalcomp && !is.null(x$B) && !is.null(x$A)) {
    As = numeric(nrow(x$sign))
    Bs = numeric(nrow(x$sign))
    for(i in 1:nrow(x$sign)) {
      As[i] = x$A[i,x$sign$index[i]]
      Bs[i] = x$B[i,x$sign$index[i]]
    }
    y = cbind(x$sign, As,Bs)
    cols = c(ncol(y)-1, ncol(y),cols)
    names(y) = c(names(x$sign),"A","B")
  }
  else y = x$sign
  
  
  #Show the number of species associated to each level
  for (k in 1:(ncolsign - 4)) {
    cat(" Number of species associated to", k, if(k==1) "group:" else "groups:", sum(rowSums(a[, 1:(ncolsign-3)]) == k), "\n")
  }
  
  cat("\n List of species associated to each combination: \n")
  
  #Show indicators for all combinations  
  for (i in 1:ncomb) {
    sel = x$sign$index == i & !is.na(x$sign$p.value) & x$sign$p.value <= alpha
    if(!is.null(minstat)) sel = sel & (abs(x$sign$stat) >= minstat)
    if(!is.null(Bt) && !is.null(x$B)) {
      for(j in 1:nrow(x$sign)) sel[j] = sel[j] && (x$B[j,x$sign$index[j]]>=Bt)
    }
    if(!is.null(At) && !is.null(x$A)) {
      for(j in 1:nrow(x$sign)) sel[j] = sel[j] && (x$A[j,x$sign$index[j]]>=At)
    }
    m = y[sel, ]
    if (nrow(m) > 0) {
      cat("\n Group", colnames(x$comb)[i], " #sps. ", nrow(m), "\n")
      m = m[order(abs(m$stat), decreasing = TRUE), cols]
      printCoefmat(m, signif.stars = TRUE, signif.legend = FALSE, 
                   digits = 4, P.values = TRUE, has.Pvalue = TRUE)
    }
  }
  Signif <- symnum(x$sign$p.value, corr = FALSE, na = FALSE, 
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))
  cat("---\nSignif. codes: ", attr(Signif, "legend"), "\n")
}
