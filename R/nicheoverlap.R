#' Metrics to compare pairs of resource niches
#' 
#' Functions \code{nicheoverlap} and \code{nichedispl} compute the overlap and centroid distance between pairs of resource distributions. In both cases resource relationships are given in the distance matrix \code{D} and the resource use data are given in data frame \code{P1} (and in some modes also \code{P2}).
#'
#' @param P1 Data frame containing the amount of usage that a set of species (in rows) make of a first set of resources (in columns)
#' @param P2 Data frame containing the amount of usage that a set of species (in rows) make of a second set of resources (in columns). Not used if \code{mode = "pairwise"}
#' @param D Object of type \code{\link{dist}} containing distance values between resources. If no distance matrix is provided (i.e. if \code{D==NULL}), the distances between resources is assumed to be maximum
#' @param q1 Vector with the availability of each resource corresponding to P1
#' @param q2 Vector with the availability of each resource corresponding to P2
#' @param mode Either \code{mode = "single"} (rows of matrices P1 and P2 are individual observations to be pooled, for example to compare the niche of two species each with its individual observations), \code{mode = "multiple"} (each row in P1 is compared to the corresponding row of P2, for example, to compare seasonal niche shifts in each species) or \code{mode = "pairwise"} (all rows in P1 are compared pairwise).
#' @param Np1 Vector with the number of observations per species from which the values in \code{P1} come (in \code{mode = "multiple"} or \code{mode = "pairwise"}).
#' @param Np2 Vector with the number of observations per species from which the values in \code{P2} come (only in \code{mode = "multiple"}).
#' @param Nq1 The number of observations from which the values in \code{q1} come.
#' @param Nq2 The number of observations from which the values in \code{q2} come.
#' @param nboot Number of boostrap samples used to compute bias-corrected percentile confidence intervals.
#' @param alpha Used to set the confidence level (i.e. \code{alpha = 0.05} means 95 percent confidence interval).
#'
#' @details
#' The method is described in De Caceres et al. (2011). If the distance matrix is not specified (i.e. if \code{D=NULL}) the function assumes that all resources are at a maximum distance (\code{d=1}). If the resource availability vector \code{q1} (and \code{q2} if supplied) is specified, then the values in \code{P1} (and \code{P2} if supplied) are taken as assessments of resource use and the species preference is calculated taking into account resource availability. Otherwise, resource use is equated to resource preference (i.e. all resources are considered equally available). The functions can compute bootstrap confidence intervals following the bias-corrected percentile method (Manly 2007). If \code{mode = "multiple"} and \code{Np1} and \code{Np2} are not null, bootstrap samples for a given niche are generated assuming a multinomial distribution with the proportions calculated from the corresponding row values in \code{P1} (resp. \code{P2}), and the number of observations comes from the corresponding element in \code{Np1} (resp. \code{Np2}). Similarly, if \code{mode = "pairwise"} and \code{Np1} is not null, bootstrap samples for each niche are generated assuming a multinomial distribution with the proportions calculated from the corresponding row values in \code{P1}, and the number of observations comes from the corresponding element in \code{Np1}. Finally, if \code{mode = "single"} then the bootstrapped units are the rows of matrices \code{P1} and  \code{P2}. In both cases, if \code{Nq1} (and \code{Nq2}) is indicated, the availability of resources is also bootstrapped. The bias-corrected percentile method is described for overlap niche measures in Mueller and Altenberg (1985).
#' 
#' @return 
#' Function \code{nicheoverlap} (resp. \code{nichedispl}) returns the overlap (resp. the distance between centroids) between the each pair of rows in \code{P1} and \code{P2}. If \code{mode = "multiple"} or \code{mode = "single"} the values are returned as a data frame. If \code{mode = "pairwise"} a matrix of values is returned instead. If bootstrap confidence intervals are asked then the functions also compute the lower and upper bounds of a confidence interval obtained following the bias-corrected percentile method. Upper and lower bounds are returned as additional columns of the data frame in \code{mode = "multiple"} or \code{mode = "single"} or as additional matrices of a list in \code{mode = "pairwise"}.
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
#' @seealso See \code{\link{nichevar}} for descriptors of single niches.
#' 
#' @export
#'
#' @name nicheoverlap
#' @examples
#' # Loads example data
#' data(birds)
#' 
#' # The overlap and displacement metrics using distances among 
#' # resources and assuming equal availability of resources
#' nicheoverlap(birdsbreed, birdswinter, D = resourceD, mode="multiple") 
#' nichedispl(birdsbreed, birdswinter, D = resourceD, mode="multiple") 
#' 
#' # The overlap and displacement metrics using distances among resources
#' # and computes 95 percent confidence intervals
#' nicheoverlap(birdsbreed, birdswinter, D = resourceD, mode="multiple", 
#'              Np1 = rowSums(birdsbreed), Np2 = rowSums(birdswinter), Nq1 = 100, Nq2 = 100) 
#' nichedispl(birdsbreed, birdswinter, D = resourceD, mode="multiple", 
#'            Np1 = rowSums(birdsbreed), Np2 = rowSums(birdswinter), Nq1 = 100, Nq2 = 100) 
#' 
#' # Same computations with different resource availability
#' q = c(0.18, 0.24, 0.22, 0.21, 0.15)
#' nicheoverlap(birdsbreed, birdswinter, D = resourceD, 
#'              q1 = q, q2 = q, mode="multiple")
#' nichedispl(birdsbreed, birdswinter, D = resourceD, 
#'            q1 = q, q2 = q, mode="multiple")
#' nicheoverlap(birdsbreed, birdswinter, D = resourceD, 
#'              q1 = q, q2 = q, mode="multiple", 
#'              Np1 = rowSums(birdsbreed), Np2 = rowSums(birdswinter), 
#'              Nq1 = 100, Nq2 = 100)
#' nichedispl(birdsbreed, birdswinter, D = resourceD, 
#'            q1 = q, q2 = q, mode="multiple", 
#'            Np1 = rowSums(birdsbreed), Np2 = rowSums(birdswinter), 
#'            Nq1 = 100, Nq2 = 100)
#' 
#' # The overlap metrics using distances among rows of 'birdsbreed'
#' nicheoverlap(birdsbreed, D = resourceD, mode="pairwise") 
nicheoverlap <- function (P1, P2 = NULL, D = NULL, q1 = NULL, q2 = NULL, mode = "multiple", Np1 = NULL, 
            Np2 = NULL, Nq1 = NULL, Nq2 = NULL, nboot = 1000, alpha = 0.05) {
  MODES <- c("single", "multiple", "pairwise")
  mode <- match.arg(mode, MODES)
  if(mode =="multiple" || mode =="single") {
    if(is.null(P2)) stop("P2 cannot be null in mode 'multiple' or 'single'")
    if (!inherits(P1, "data.frame") || !inherits(P2, "data.frame")) stop("P1 and P2 should be dataframes")
    if (mode == "multiple" && (nrow(P1) != nrow(P2))) stop("Resource use dataframes do not have the same number of rows")
    if (ncol(P1) != ncol(P2)) stop("Resource use dataframes do not have the same number of columns (resources)")  
    if (!is.null(Np1) && !is.null(Np2)) {
      if (!inherits(Np1, "numeric") || !inherits(Np2, "numeric")) stop("Np1 and Np2 should be numeric vectors")
      Np1 = as.integer(Np1)
      Np2 = as.integer(Np2)
    }
    if((!is.null(Nq1) && is.null(Nq2)) || (is.null(Nq1) && !is.null(Nq2))) stop("Nq1 and Nq2 should be both either NULL or contain numeric values")
    if (!is.null(Nq1) && !is.null(Nq2)) {
      if (!inherits(Nq1, "numeric") || !inherits(Nq2, "numeric")) stop("Nq1 and Nq2 should be numeric")
      Nq1 = as.integer(Nq1)
      Nq2 = as.integer(Nq2)
    }
    if (!is.null(D)) {
      if (!inherits(D, "dist")) 
        stop("Object of class 'dist' expected for distance")
      D <- as.matrix(D)
      if (ncol(P1) != nrow(D)) stop("The number of columns in P1 must be equal to the number of items in D")
      if (ncol(P2) != nrow(D)) stop("The number of columns in P2 must be equal to the number of items in D")
      D <- as.dist(D)
    }
  } else if(mode=="pairwise") {
    if (!inherits(P1, "data.frame")) stop("P1 should be a dataframe")
    if (!is.null(D)) {
      if (!inherits(D, "dist")) stop("Object of class 'dist' expected for distance")
      D <- as.matrix(D)
      if (ncol(P1) != nrow(D)) stop("The number of columns in P1 must be equal to the number of items in D")
      D <- as.dist(D)
    }    
  }
  #Check 'q1' the availability of resources in the first 'season'
  if (!is.null(q1)) {
    if (length(q1) != ncol(P1)) stop("The number of items in q1 must be equal to the number of columns in P1 and P2")
    q1 = q1/sum(q1)
  } else {
    q1 = rep(1/ncol(P1), ncol(P1))
  }
  #Check 'q2' the availability of resources in the second 'season'
  if (!is.null(q2)) {
    if (length(q2) != ncol(P2)) stop("The number of items in q1 must be equal to the number of columns in P1 and P2")
    q2 = q2/sum(q2)
  } else {
    q2 = rep(1/ncol(P2), ncol(P2))
  }
  
  #If no distance matrix was supplied, generate one (equidistant resources)
  if (is.null(D)) D <- as.dist((matrix(1, ncol(P1), ncol(P1)) - diag(rep(1, ncol(P1)))))
  
  #Internal functions
  nichevar1 <- function(f, D) {
      if (is.na(sum(f))) v <- NA
      else if (sum(f) < 1e-16) v <- 0
      else v <- (f %*% (as.matrix(D)^2) %*% f)/(2 * (sum(f)^2))
      return(v)
  }
  getF <- function(p, q = NULL) {
    if (!is.null(q)) { 
      a = p/q
      return(a/sum(a))
    } else {
      return(p/sum(p))
    }
  }
  overlap1 <- function(f1, f2, D) {
    if (is.na(sum(f1)) || is.na(sum(f2))) o <- NA
    else if (sum(f1) < 1e-16 || sum(f2) < 1e-16) o <- 0
    else {
      o <- (1 - ((f1 %*% (as.matrix(D)^2) %*% f2)/(sum(f2) * sum(f1))))/sqrt((1 - 2 * nichevar1(f1, D)) * (1 - 2 * nichevar1(f2, D)))
    }
    return(o)
  }
  
    
  #Calculate overlap between each row of P1 and the corresponding row in P2
  if (mode == "multiple") {
    if (!is.null(Np1) && !is.null(Np2)) nc = 3 #If we have to calculate confidence intervals
    else nc = 1
    O <- as.data.frame(matrix(0, nrow = nrow(P1), ncol = nc))
    for (i in 1:nrow(P1)) rownames(O)[i] <- paste(row.names(P1)[i], "vs", row.names(P2)[i])
    for (i in 1:nrow(P1)) {
      pi1 = as.numeric(P1[i, ])
      pi2 = as.numeric(P2[i, ])
      O[i, 1] <- overlap1(getF(pi1, q1), getF(pi2, q2), D)
      if (!is.null(Np1) && !is.null(Np2)) {
        BO = vector("numeric", length = nboot)
        O[i, 2] = NA
        O[i, 3] = NA
        if (!is.na(sum(pi1)) && !is.na(sum(pi2))) {
          bsamp1 = rmultinom(nboot, Np1[i], getF(pi1))
          bsamp2 = rmultinom(nboot, Np2[i], getF(pi2))
          if (!is.null(Nq1) && !is.null(Nq2)) {
            qsamp1 = rmultinom(nboot, Nq1, q1)
            qsamp2 = rmultinom(nboot, Nq2, q2)
          }
          for (b in 1:nboot) {
            if (!is.null(Nq1) && !is.null(Nq2)) BO[b] = overlap1(getF(bsamp1[, b], qsamp1[, b]), getF(bsamp2[, b], qsamp2[, b]), D)
            else BO[b] = overlap1(getF(bsamp1[, b], q1), getF(bsamp2[, b], q2), D)
          }
          #Some NA may appear because of zeroes in qsamp
          BO = BO[!is.na(BO)]
          #Compute Bias-corrected percentile method (Manly 2007: pp52-56)
          z0 = qnorm(sum(BO < O[i, 1])/length(BO))
          lj = floor(length(BO) * pnorm(2 * z0 + qnorm(alpha/2)))
          uj = floor(length(BO) * pnorm(2 * z0 + qnorm(1 - (alpha/2))))
          if (lj > 0 && uj > 0 && lj != uj) {
            sbo = sort(BO)
            O[i, 2] = sbo[lj]
            O[i, 3] = sbo[uj]
          }
        }
      }
    }
    if (nc == 1) names(O) <- "O"
    else names(O) <- c("O", "LC", "UC")
    return(O)
  }
  
  #P1 are different resource use assessments of a single entity (the same for P2)
  if (mode == "single") {
    O <- as.data.frame(matrix(NA, nrow = 1, ncol = 3))
    rownames(O) <- "Overlap"
    O[1, 1] <- overlap1(getF(colSums(P1, na.rm=TRUE), q1), getF(colSums(P2, na.rm=TRUE), q2), D)
    BO = vector("numeric", length = nboot)
    if (!is.null(Nq1) & !is.null(Nq2)) {
      qsamp1 = rmultinom(nboot, Nq1, q1)
      qsamp2 = rmultinom(nboot, Nq2, q2)
    }
    for (b in 1:nboot) {
      p1samp = colSums(P1[sample(1:nrow(P1), replace = TRUE), ], na.rm=TRUE)
      p2samp = colSums(P2[sample(1:nrow(P2), replace = TRUE), ], na.rm=TRUE)
      if (!is.null(Nq1) & !is.null(Nq2)) BO[b] = overlap1(getF(p1samp, qsamp1[, b]), getF(p2samp, qsamp2[, b]), D)
      else BO[b] = overlap1(getF(p1samp, q1), getF(p2samp, q2), D)
    }
    BO = BO[!is.na(BO)]
    z0 = qnorm(sum(BO < O[1, 1])/length(BO))
    lj = floor(length(BO) * pnorm(2 * z0 + qnorm(alpha/2)))
    uj = floor(length(BO) * pnorm(2 * z0 + qnorm(1 - (alpha/2))))
    if (lj > 0 && uj > 0 && lj != uj) {
      sbo = sort(BO)
      O[1, 2] = sbo[lj]
      O[1, 3] = sbo[uj]
    }
    names(O) <- c("O", "LC", "UC")
    return(O)
  }
  
  #Only P1 is used. Calculate overlap between pairs of rows in P1. No confidence intervals are calculated
  if(mode=="pairwise") {
    O <- matrix(1, nrow = nrow(P1), ncol = nrow(P1))
    rownames(O)<-rownames(P1)
    colnames(O)<-rownames(P1)
    if(!is.null(Np1)) {
      LC <- O
      UC <- O
    }
    for (i in 1:(nrow(P1)-1)) {
      for (j in (i+1):nrow(P1)) {
        pi = as.numeric(P1[i, ])
        pj = as.numeric(P1[j, ])
        O[i,j] <- overlap1(getF(pi, q1), getF(pj, q1), D)
        O[j,i] <- O[i,j]
        if (!is.null(Np1)) {
          BO = vector("numeric", length = nboot)
          if (!is.na(sum(pi)) && !is.na(sum(pj))) {
            bsampi = rmultinom(nboot, Np1[i], getF(pi))
            bsampj = rmultinom(nboot, Np1[j], getF(pj))
            if (!is.null(Nq1)) qsamp1 = rmultinom(nboot, Nq1, q1)
            for (b in 1:nboot) {
              if (!is.null(Nq1)) BO[b] = overlap1(getF(bsampi[, b], qsamp1[, b]), getF(bsampj[, b], qsamp1[, b]), D)
              else BO[b] = overlap1(getF(bsampi[, b], q1), getF(bsampj[, b], q1), D)
            }
            BO = BO[!is.na(BO)]
            z0 = qnorm(sum(BO < O[i,j])/length(BO))
            lj = floor(length(BO) * pnorm(2 * z0 + qnorm(alpha/2)))
            uj = floor(length(BO) * pnorm(2 * z0 + qnorm(1 - (alpha/2))))
            if (lj > 0 && uj > 0 && lj != uj) {
              sbo = sort(BO)
              LC[i, j] = LC[j, i] =sbo[lj]
              UC[i, j] = UC[j, i] =sbo[uj]
            }
          }
        }
      }
    }
    if(!is.null(Np1)) {
      return(list(O=O, LC=LC, UC = UC))
    }
    else return(O)
  }
}

