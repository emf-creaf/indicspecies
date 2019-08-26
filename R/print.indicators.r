print.indicators <- function (x, At = 0, Bt = 0, sqrtIVt = 0, alpha = 1.0, selection = NULL, confint=FALSE,...) {
	if(is.null(selection)) selection = rep(TRUE, nrow(x$C))
	if(length(dim(x$A))==2) {
		A = x$A[selection,]
		B = x$B[selection,]
		sqrtIV = x$sqrtIV[selection,]
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
