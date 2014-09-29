# update the detVx by using the result of "split2cls()"
# continue: no update
# org.detVx: original detVx vector
# divided.detVx: divided detVx vector having 2 elements.
# Author: tunenori
###############################################################################

updtdetVx <- function(continue, org.detVx, k1, k2, divided.detVx){
	if (!is.vector(org.detVx))
		return(divided.detVx)
	if (!continue)
		return(org.detVx)
	if (k1 == k2)
		stop("updtdetVx() : k1 and k2 should differ.")
	
	z <- NULL
	for (i in 1:max(k2, length(org.detVx))){
		if (i == k1)
			z <- c(z, divided.detVx[1])
		else if (i == k2)
			z <- c(z, divided.detVx[2])
		else
			z <- c(z, org.detVx[i])
	}
	z
}

