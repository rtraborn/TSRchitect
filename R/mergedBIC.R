# return BICs for Two-merged clusters model and devided clusters model
# k1/k2: marged cluster ID
# Author: tunenori
###############################################################################

mergedBIC <- function(x, xcl, k1, k2, q, ignore.covar, pr.proc){
	# sample size
	# check for input data
	n1 <- xcl$size[k1]
	n2 <- xcl$size[k2]
	if (n1 == 0 || n2 == 0){
		# already had been merged
		cat(paste("already had been merged\n"))
		ret <- F
		return( list (ret = ret))
	}
	if (is.null(xcl$lnL[k1]) || is.null(xcl$lnL[k2])){
		# lnL may be null because of few data
		cat(paste("lnL may be null because of few data\n"))
		ret <- F
		return( list (ret = ret))
	}
	
	# divided clusters model
	lnL1 = xcl$lnL[k1]
	lnL2 = xcl$lnL[k2]
	ctrextrt <- rbind(xcl$centers[k1,], xcl$centers[k2,])
	beta <- dist(ctrextrt) / (sqrt(xcl$detVx[k1] + xcl$detVx[k2]))
	if (pr.proc) cat(paste("beta=", round (beta, digit=2), "\n"))
	
	# if (beta > 10){
	# 	# 2 clusters far apart
	# 	ret <- F
	# 	return( list (ret = ret))
	# }
	
	alpha <- 0.5 / as.numeric(pnorm(beta))
	bicdiv <- -2 * lnL1  -2 * lnL2 + 2 * q * log(n1 + n2) - 2 * (n1 + n2) * log(alpha)
	# bicdiv <- -2 * lnL1 -2 * lnL2 + 2 * q - 2 * (n1 + n2) * log(alpha) # AIC
	
	# extract 2 clusters data
	y.ok1 <- lapply(xcl$cluster, "==", k1) # 1st sub-cluster or not
	y.ok2 <- lapply(xcl$cluster, "==", k2) # 2nd sub-cluster or not
	
	# extract sub data
	p = ncol(x)
	clj1 <- matrix(x[as.logical(y.ok1)], ncol=p)
	clj2 <- matrix(x[as.logical(y.ok2)], ncol=p)
	xmgd <- rbind(clj1, clj2)
	
	# merged cluster center
	ctrmgd <- (n1 * xcl$centers[k1,] + n2 * xcl$centers[k2,]) / (n1 + n2)
	lnLmgd <- lnL(xmgd, ctrmgd, ignore.covar)
	bicmgd <- -2 * lnLmgd$lnL + q * log(nrow(xmgd)) # BIC
	# bicmgd <- -2 * lnLmgd$lnL + q  # AIC
	
	ret <- T
	list (ret = ret, ctrmgd = ctrmgd, lnLmgd = lnLmgd$lnL, detVxmgd = lnLmgd$detVx, bicmgd = bicmgd, bicdiv = bicdiv)
}
