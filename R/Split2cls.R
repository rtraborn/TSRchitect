# split 2 clusters if we would prefer it based on BIC
# q: a number of parameters
# bic.prior: BIC which x is given; if bic.prior=NULL then we calculate
# lnL.prior: lnL which x is given; if bic.prior=NULL then we calculate
# detVx.prior: detVx which x is given; if bic.prior=NULL then we calculate
# Author: tunenori
###############################################################################

split2cls <- function(x, centers, q, bic.prior, lnL.prior, detVx.prior, iter.max, ignore.covar){
	if (is.null(bic.prior)){
		pb <- priorBIC(x, centers, q, ignore.covar)
		bic.prior <- pb$bic
		lnL.prior <- pb$lnL
		detVx.prior <- pb$detVx
	}
	bic.post <- postBICs(x, centers, q, iter.max, ignore.covar)
	
	subcluster <- bic.post$clsub$cluster
	
	# compare whether if we should split
	if (is.na(bic.post$bic[3])){
		# BIC may has NA because of few data 
		continue <- FALSE
	}else if (bic.post$bic[3] < bic.prior){
		# splitting ...
		# replace the cluster number to cl$cluster
		continue <- TRUE
	}else{
		# not splitting...
		# return "subcluster" stored k1 
		continue <- FALSE
	}
	# note that "subcluster" gives 1 or 2 
	list(continue = continue, subcluster = subcluster, 
			bic.prior = bic.prior, bic.post = bic.post$bic,
			lnL.prior = lnL.prior, lnL.post = bic.post$lnL,
			detVx.prior = detVx.prior, detVx.post = bic.post$detVx,
			centers = bic.post$clsub$centers,
			clj1 = bic.post$clj1, clj2 = bic.post$clj2)
}
