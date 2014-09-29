# return BICs (two posterior BICs)
# Author: tunenori
###############################################################################

postBICs <- function(x, centers, q, iter.max, ignore.covar){
	#
	# split to 2 clusters
	clsub <- kmeans_modv2(x, 2, iter.max, algorithm = "Lloyd")
	y.ok1 <- lapply(clsub$cluster, "==", 1) # 1st sub-cluster or not
	y.ok2 <- lapply(clsub$cluster, "==", 2) # 2nd sub-cluster or not
	# extract sub data
	p <- ncol(x)
	clj1 <- matrix(x[as.logical(y.ok1)], ncol=p)
	clj2 <- matrix(x[as.logical(y.ok2)], ncol=p)
	# ratio for pdf.
	r1 <- clsub$size[1] / sum(clsub$size)	# [0,1]
	r2 <- 1 - r1 	# [0,1]
	# two later BICs
	# print(clsub$centers[1,])	# for debug
	# print(apply(clj1,2,mean))	# for debug
	# print(sqrt(apply(clj1,2,var)))	# for debug
	# print(r1)	# for debug
	lnL1 <-  lnL(clj1, clsub$centers[1,], ignore.covar)
	# print(clsub$centers[2,])	# for debug
	# print(apply(clj2,2,mean))	# for debug
	# print(sqrt(apply(clj2,2,var)))	# for debug
	# print(r2)	# for debug
	lnL2 <-  lnL(clj2, clsub$centers[2,], ignore.covar)
	n1 <- nrow(clj1)
	n2 <- nrow(clj2)
	# normalizing factor; dist() is in library(mva)
	if (is.na(lnL1$detVx) || is.na(lnL2$detVx))
		beta <- 0
	else
		beta <- dist(clsub$center) / (sqrt(lnL1$detVx + lnL2$detVx))
	alpha <- 0.5 / pnorm(beta)
	BIC1 <- -2 * lnL1$lnL +q * log(n1)
	BIC2 <- -2 * lnL2$lnL +q * log(n2) 
	# BIC1 <- -2 * lnL1$lnL +q # AIC
	# BIC2 <- -2 * lnL2$lnL +q # AIC
	
	# cat (paste("alpha =",alpha,"\n"))	# for debug
	# cat (paste("beta =",beta,"\n"))	# for debug
	
	# BIC is not (BIC1 + BIC2)
	BIC <- -2 * lnL1$lnL  -2 * lnL2$lnL + 2 * q * log(n1 + n2) - 2 * (n1 + n2) * log(alpha)
	# BIC <- -2 * lnL1$lnL  -2 * lnL2$lnL + 2 * q  - 2 * (n1 + n2) * log(alpha) # AIC
	list(bic = c(BIC1, BIC2, BIC), 
			lnL = c(lnL1$lnL, lnL2$lnL),
			detVx = c(lnL1$detVx, lnL2$detVx),
			clsub = clsub, clj1 = clj1, clj2 = clj2)
}


