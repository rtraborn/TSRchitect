# update the lnL by using the result of "split2cls()"
# continue: no update
# org.lnL: original lnL vector
# divided.lnL: divided lnL vector having 2 elements.
# Author: tunenori
###############################################################################

updtlnL <- function(continue, org.lnL, k1, k2, divided.lnL){
	if (!is.vector(org.lnL))
		return(divided.lnL)
	if (!continue)
		return(org.lnL)
	if (k1 == k2)
		stop("updtlnL() : k1 and k2 should differ.")
	
	z <- NULL
	for (i in 1:max(k2, length(org.lnL))){
		if (i == k1)
			z <- c(z, divided.lnL[1])
		else if (i == k2)
			z <- c(z, divided.lnL[2])
		else
			z <- c(z, org.lnL[i])
	}
	z
}
