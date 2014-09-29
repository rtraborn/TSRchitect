# update the cluster centers by using the result of "split2cls()"
# continue: no update
# org.centers: original centers matrix
# divided.centers: divided centers matrix; it has 2 rows.
# Author: tunenori
###############################################################################

updtCenters <- function(continue, org.centers, k1, k2, divided.centers){
	if (!is.matrix(org.centers)) 
		return(divided.centers)
	if (!continue)
		return(org.centers)
	if (k1 == k2)
		stop("updtCenters() : k1 and k2 should differ.")
	
	z <- NULL
	for (i in 1:max(k2, nrow(org.centers))){
		if (i == k1)
			z <- rbind(z, divided.centers[1,])
		else if (i == k2)
			z <- rbind(z, divided.centers[2,])
		else
			z <- rbind(z, org.centers[i,])
	}
	z
}

