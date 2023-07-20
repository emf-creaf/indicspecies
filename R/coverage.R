#' Coverage of a set of indicators
#' 
#' Function \code{coverage} calculates the proportion of sites of the target site group where one or another indicator (a species or a species combination) is found. Parameters are used to select valid indicators before calculating coverage. Function \code{plotcoverage} plots the coverage against the positive predictive value threshold used to select indicators.
#'  
#' @param x An object of class '\code{\link{indicators}}' or a community data table with sites in rows and species in columns (in this second case, an object of class '\code{\link{multipatt}}' must be supplied for \code{y}).
#' @param y An object of class '\code{\link{multipatt}}'.
#' @param selection A logical vector restricting the set of indicators used to calculate the coverage.
#' @param minstat Minimum value of the statistic for selecting indicators.
#' @param At Minimum value of positive predictive value (A) for selecting indicators.
#' @param Bt Minimum value for sensitivity (B) for selecting indicators.
#' @param type Specifies how to select indicators: (1) using the value of the statistic (\code{type = "stat"}); 
#' (2) the lower bound of its confidence interval (\code{type = "lowerCI"}); 
#' or (3) the upper bound of its confidence interval (\code{type = "upperCI"}). 
#' This parameter makes sense when the function is called using objects of class '\code{indicators}' 
#' and bootstrap confidence intervals are available for this object. Otherwise \code{type} has no 
#' effect and the value of the statistic is used for selection. In function \code{coverage}, 
#' the value of \code{type} applies to selection using \code{minstat}, \code{At} and \code{Bt}. 
#' In function \code{plotcoverage}, the value of \code{type} applies to selection using \code{At}.
#' @param alpha Significance level for selecting indicators.
#' @param by Rate of increase in the predictive value threshold (At).
#' @param max.order The maximum number of species conforming species combinations (allows examining the effects of increasing the order of combinations).This parameter is only used when the function is called using objects of class '\code{indicators}'.
#' @param group Either an integer or a character string indicating the site group or site group combination for which plot is desired. This parameter is only used when the function is called using objects of class '\code{multipatt}'.
#' @param add Flag indicating whether the function should draw on the existing plot.
#' @param xlab Label for the x-axis.
#' @param ... Additional plotting parameters that are passed to the \code{plot} function.
#'
#' @details
#' The \code{coverage} of a set of indicators was defined in De \enc{Cáceres}{Caceres} et al. (2012) as the proportion of sites in a given site group where one or several indicators are found. This value allows assessing how often the site group will be able to be determined. If all indicators of a site group are rare, then the indication system will not be useful, regardless of how much restricted to the site group the indicators are. The coverage value is a generalization of quantity B of IndVal, that applies to a group of indicators instead of a single one. Function \code{plotcoverage} plots the coverage against the positive predictive value threshold (At) used to select indicators, as in De \enc{Cáceres}{Caceres} et al. (2012). Functions \code{coverage} and \code{plotcoverage} can be executed using either an object of class '\code{indicators}', or an object of class '\code{multipatt}'. However, the parameters that apply to each case are slightly different. When using \code{coverage} and \code{plotcoverage} on objects of class '\code{multipatt}' one is expected to calculate the coverage for those indicators that are significant (see \code{alpha} parameter), although other constraints to select valid indicators can be used. When using \code{coverage} and \code{plotcoverage} on objects of class '\code{indicators}' one is expected to calculate the coverage for indicators that have values of A larger than a specified threshold (see \code{At} parameter). In this latter case, it may be advisable to use \code{stat="lowerCI"}, so that indicators with broad confidence intervals are not included in the selection.
#' 
#' @return
#' When used with an object of class '\code{indicators}', function \code{coverage} returns the proportion of sites of the target site group where one or another indicator (species combination) is found. When used with an object of class '\code{indicators}', function \code{coverage} returns a vector containing the coverage value for each site group or site group combination.
#' 
#' @references 
#' De \enc{Cáceres}{Caceres}, M., Legendre, P., Wiser, S.K. and Brotons, L. 2012. Using species combinations in indicator analyses. Methods in Ecology and Evolution 3(6): 973-982.
#' 
#' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, EMF-CREAF
#' 
#' @seealso \code{\link{indicators}}, \code{\link{multipatt}}, \code{\link{pruneindicators}}
#' 
#' @export
#' 
#' @name coverage
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
#' ## Determine the coverage of the selected set of indicators
#' coverage(sc)
#' 
#' ## Plot the coverage against the threshold At
#' plotcoverage(sc)
#' plotcoverage(sc, max.order=2, add=TRUE, lty=2)
#' 
#' ## Runs the combination analysis using IndVal.g as statistic
#' wetpt <- multipatt(wetland, wetkm$cluster, control = how(nperm=999))
#'  
#' ## Determines the coverage for each site group combination
#' coverage(wetland, wetpt, alpha = 0.05)
#' 
coverage <- function (x, y=NULL, selection=NULL, minstat=NULL, At=NULL, Bt=NULL, type="stat", alpha=NULL) {
    match.arg(type,c("lowerCI","upperCI","stat"))
    if(inherits(x,"indicators")) {
      speciescomb = x
      if(is.null(selection)) selection = rep(TRUE, nrow(speciescomb$C))
      if(!is.null(At)) {
        if(is.data.frame(speciescomb$A)) {
          if(type=="stat") selection = selection & (speciescomb$A$stat>=At)
          else if(type=="lowerCI") selection = selection & (speciescomb$A$lowerCI>=At)
          else if(type=="upperCI") selection = selection & (speciescomb$A$upperCI>=At)
        }
        else selection = selection & (speciescomb$A>=At)
      }
      if(!is.null(Bt)) {
        if(is.data.frame(speciescomb$B)) {
          if(type=="stat") selection = selection & (speciescomb$B$stat>=Bt)
          else if(type=="lowerCI") selection = selection & (speciescomb$B$lowerCI>=Bt)
          else if(type=="upperCI") selection = selection & (speciescomb$B$upperCI>=Bt)
        }
        else selection = selection & (speciescomb$B>=Bt)
      }
      if(!is.null(minstat)) {
        if(is.data.frame(speciescomb$sqrtIV)) {
          if(type=="stat") selection = selection & (speciescomb$sqrtIV$stat>=minstat)
          else if(type=="lowerCI") selection = selection & (speciescomb$sqrtIV$lowerCI>=minstat)
          else if(type=="upperCI") selection = selection & (speciescomb$sqrtIV$upperCI>=minstat)
        }
        else selection = selection & (speciescomb$sqrtIV>=minstat)
      }
      if(!is.null(alpha)) {
        selection = selection & (speciescomb$p.value<=alpha)
      }
      
      if(length(dim(speciescomb$C))==2) c = speciescomb$C[selection,]
      else c = speciescomb$c[selection]
      if(sum(selection)==0) return(0)
      group.vec = speciescomb$group.vec
      if(length(dim(speciescomb$XC))==2) xc = speciescomb$XC[, selection]
      else xc = speciescomb$XC[selection]
      if(sum(selection)>1) {
        ccx = rep(FALSE, nrow(xc))    
        for(rc in 1:nrow(c)) {
          ccx = ccx | xc[,rc]>0  		
        }
        return(sum(ccx & group.vec) / sum(group.vec))
      } else {
        return(sum(xc>0 & group.vec) / sum(group.vec))
      }      
    }
    else if(inherits(x,"data.frame")) {
      if(is.null(y)) stop("You must supply a multipatt object in 'y'")
      if(!inherits(y,"multipatt")) stop("Wrong class for 'y'. Should be an object of class 'multipatt'")
      mp = y
      ncomb = ncol(mp$comb)
      cov = numeric(ncomb)
      for(c in 1:ncomb) {
        ind = mp$sign$index
        sel = ind == c 
        selgroup = as.logical(mp$comb[,c])
        if(!is.null(selection)) sel = sel & selection
        if(!is.null(At) && (mp$func=="IndVal" || mp$func=="IndVal.g")) { #For each indicator, check A value
          selind = which(sel)
          for(i in selind) sel[i] = (mp$A[i,ind[i]]>=At) 
        }
        if(!is.null(Bt) && (mp$func=="IndVal" || mp$func=="IndVal.g")) { #For each indicator, check B value
          selind = which(sel)
          for(i in selind) sel[i] = (mp$B[i,ind[i]]>=Bt) 
        }
        if(!is.null(minstat)) { #For each indicator, check minstat
          sel = sel & (mp$sign$stat>=minstat)
        }
        if(!is.null(alpha)) { #For each indicator, check p-value
          sel = sel & (mp$sign$p.value<=alpha)
          sel[is.na(sel)] = FALSE #There may be NA values for the combination of all sites
        }
        if(sum(sel)==1) cov[c] = sum(x[selgroup,sel]>0)/sum(selgroup)
        else {
          cov[c] = sum(rowSums(x[selgroup,sel])>0)/sum(selgroup)
        }
      }
      names(cov) = colnames(mp$comb)
      return(cov)
    }
}