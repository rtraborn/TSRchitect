# update the cluster number by using the result of "split2cls()"
# continue: no splitting
# v: cluster numbers vector for initial cluster.
# k1: cluster numbers should be updated; "k1" becomes "k1" and "k2"
# xsub: sub-cluster numbers vector of "v" whose value is "k1";
# given "xsub" have 1 or 2.
# Author: tunenori
###############################################################################

updtCrusterNum <- function(continue, v, k1, k2, xsub){
	if (!is.vector(v)) 
		return(xsub)
	if (!continue)
		return(v)
	if (k1 == k2)
		stop("updtCrusterNum() : k1 and k2 should differ.")
	
	# below is same algorithm; explicit array operation is slow in R.
	# j <- 1
	# for (i in 1:length(v)){
	#	if (v[i] == k1){
	#		if (xsub[j] == 2)
	#			v[i] <- k2
	#		j <- j + 1
	#	}
	# }
	# end of algorithm
	xsub <- replace(xsub, (xsub == 2), k2) # changed
	xsub <- replace(xsub, (xsub == 1), k1) # unchanged
	v <- replace(v, (v == k1), xsub)
}

