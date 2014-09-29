# marge the result of sub-clustering;
# cluster numbers by first kmeans should be renumbered;
# the other centers and sizes are simply added.
# cl: the result of first kmeans
# cl.sub: the result of subclustering
# ik: cluster number adopted to kmeans.
# Author: tunenori
###############################################################################

mergeResult <- function(cl, cl.sub, ik){
	cluster <- cl$cluster	# main cluster
	centers <- NULL
	size <- NULL
	lnL <- NULL
	detVx <- NULL
	
	k <- 0	# uniq cluster numbers; k should be decremental. 
	for (i in 1:ik)
		k <- k + length(cl.sub[[i]]$size)
	kk <- k
	
	for (i in ik:1){	# loop for main clusters obtained by kmeans
		xsub <- cl.sub[[i]]$cluster
		iki <- ik -i +1 
		centers <- rbind(centers, cl.sub[[iki]]$centers)
		size <- c(size, cl.sub[[iki]]$size)
		lnL <- c(lnL, cl.sub[[iki]]$lnL)
		detVx <- c(detVx, cl.sub[[iki]]$detVx)
		
		for (j in length(cl.sub[[i]]$size):1){ # loop for subclusters
			xsub <- replace(xsub, (xsub == j), k)
			k <- k -1
		}
		cluster <- replace(cluster, (cluster == i), xsub)
	}
	if (k != 0) stop("mergeResult: assertion failed (k = 0)...")
	dimnames(centers) <- list(1:kk, NULL)
	list(cluster = cluster, centers = centers, lnL = lnL, detVx = detVx, size = size)
}


