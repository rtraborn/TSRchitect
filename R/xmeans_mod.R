# non-hierarchal clustering;
# it is like kmeans but cluster number is not fixed;
# a user gives initial k.
# x: data
# ik: initial k applied kmeans; k should be sufficiently small
# iter.max: maximum of iterations applied for "kmeans"
# pr.proc: print out the processing status
# ignore.covar: ignore covariance of cluster data
# merge.cls: some clusters will be merged after iterative division

# Author: tunenori
###############################################################################

xmeans_mod <- function(x, ik = 2, iter.max = 1000, pr.proc = F, ignore.covar=T, merge.cls=F){
	if (ik < 2) 
		ik <- 2
	x <- as.matrix(x)
	p <- ncol(x) # p-dimensional multivariate
	if (ignore.covar){
		q <- 2 * p # number of parameters; mean and var for each "p"
	}else{
		q <- p * (p+3) / 2	# integer
	}
	cl<- kmeans_modv2(x,ik,iter.max,algorithm = "Lloyd")
	cl.sub <- list()
	
	for (i in 1:ik){ # for each cluster
		y.ok <- (cl$cluster == i) 	# i-th cluster or not
		yi <- matrix(na.omit(x[y.ok], ncol=p)) 	# extract i-th cluster
		zi <- yi 	# save the data for graphics
		yi.centers <- cl$centers[i,]
		zi.centers <- yi.centers
		yi.cluster <- cl$cluster[(cl$cluster == i)]
		yi.cluster <- rep(1, length(yi.cluster)) 
		# sub-cluster number should begin from 1
		
		k1 <- 1		# cluster number
		k2 <- k1 + 1
		bic.prior <- NULL
		stack <- list()	# divided and unproceeded data are stacked
		lnL0 <- lnL(yi, yi.centers, ignore.covar)
		yi.lnL <- lnL0$lnL
		yi.detVx <- lnL0$detVx
		
		repeat{
			
			# go through at least 1 time; 
			# y$subcluster exist...
			if (pr.proc)	cat (paste("k1 =", k1, ", k2 =", k2,"\n"))
			if (nrow(yi) == 1){ # sample size is 1
				break
			}
			y <- split2cls(yi, yi.centers, q, bic.prior, lnL.prior, detVx.prior, iter.max, ignore.covar)
			if (y$continue){ # splitting continues 
				yi.cluster <-
						updtCrusterNum(y$continue, yi.cluster, k1, k2, y$subcluster)
				zi.centers <-
						updtCenters(y$continue, zi.centers, k1, k2, y$centers)
				yi.lnL <-
						updtlnL(y$continue, yi.lnL, k1, k2, y$lnL.post)
				yi.detVx <-
						updtdetVx(y$continue, yi.detVx, k1, k2, y$detVx.post)
			}
			
			if (pr.proc) print(y$subcluster)
			if (pr.proc){ print(y$bic.prior)
				print(y$bic.post)
				# print(y$lnL.prior)	# for debug
				# print(y$lnL.post)	# for debug
				print(y$continue) }
			# cat("zi.centers=\n")	# for debug
			# print(zi.centers)	# for debug
			if (!y$continue){	# no-node
				if ((nstack <- length(stack))){ # there are stacked data
					# extract the stacked data
					if (pr.proc)
						cat(paste("extract the stacked data (", nstack, ")...\n"))
					yi <- stack[[nstack]]$data
					yi.centers <- stack[[nstack]]$centers
					bic.prior <- stack[[nstack]]$bic
					lnL.prior <- stack[[nstack]]$lnL
					detVx.prior <- stack[[nstack]]$detVx
					k1 <- stack[[nstack]]$cls
					k2 <- k2 # unchanged
					# delete the data set
					if (nstack > 1){
						stack <- stack[1:(nstack-1)]
					}else{
						stack <- list() # no stacked data
					}
					next;
				}
				# no node and no stack
				if (pr.proc)	cat ("no node and no stack...\n")
				break;
			}
			# splitting continues...
			y1 <- y$clj1	# data
			y2 <- y$clj2
			yi.ctr1 <- y$centers[1,]	# centers
			yi.ctr2 <- y$centers[2,]
			bic.prior1 <- y$bic.post[1]	# bic
			bic.prior2 <- y$bic.post[2]
			lnL.prior1 <- y$lnL.post[1]	# lnL
			lnL.prior2 <- y$lnL.post[2]
			detVx.prior1 <- y$detVx.post[1]	# detVx
			detVx.prior2 <- y$detVx.post[2]
			
			# one-hand repeats recursively...
			yi <- y1
			yi.centers <- yi.ctr1
			bic.prior <- bic.prior1
			lnL.prior <- lnL.prior1
			detVx.prior <- detVx.prior1
			# other-hand is stacked... 
			if (pr.proc)	cat ("stacking ...\n")
			stack <- c(stack,
					list(list(data=y2, centers=yi.ctr2,
									bic=bic.prior2, lnL=lnL.prior2, detVx=detVx.prior2, cls=k2)))
			# inclement the cluster number 
			k2 <- k2 + 1
			
		} # end of repeat
		
		# splitting done ...
		if (pr.proc){
			cat ("splitting done...\n")
			cat (paste("main cluster =",i,"*******\n"))
		}
		cl.sub <- c(cl.sub, list(list(cluster = yi.cluster,
								centers = zi.centers, lnL = yi.lnL, detVx = yi.detVx,
								size = tabulate(yi.cluster))))
		if (pr.proc){
			print(cl.sub[[i]])
			plot(zi, col=yi.cluster)
			if (is.vector(zi.centers))
				points(zi.centers[1], zi.centers[2], pch=8)
			else # array
				points(zi.centers,col=1:(length(zi.centers)/p),pch=8)
		}
	}
	if (pr.proc)	print(cl.sub)
	xcl <- mergeResult(cl, cl.sub, ik)
	
	if (merge.cls == F) {
		return(list(cluster = xcl$cluster, centers = xcl$centers, size = xcl$size))
	}
	
	# merge after progressive dividing
	#
	if (pr.proc) cat("merging after progressive dividing ...\n")
	
	k <- length(xcl$size)	# final cluster number
	if (k <= 2){	# minimum cluster number should be 2
		if (pr.proc) cat("merging skipped ...\n")
		return(list(cluster = xcl$cluster, centers = xcl$centers, size = xcl$size))
	}
	if (pr.proc){
		cat("xcl$detVx=")
		print(xcl$detVx)
		cat("xcl$size=")
		print(xcl$size)
	}
	
	klist <- sort.list(xcl$size) # "small" to "large" order of xcl$detVx list
	if (pr.proc) print(klist)
	for (i in 1:(k-1)){ 
		for (j in (i+1):k){ 
			k1 = klist[i]
			k2 = klist[j]
			if (pr.proc) cat(paste("inspecting the clusters", k1,"and", k2,"\n"))
			
			z <- mergedBIC(x, xcl, k1, k2, q, ignore.covar, pr.proc)
			if (z$ret == F){
				# k1 or k2 has been merged.
				# skip this loop
				if (pr.proc) cat("skipping... k1=", k1, "k2=", k2,"\n")
				next
			}
			if (z$bicdiv > z$bicmgd){
				# we prefer merged model.
				# replace larger cls. number to smaller cls. number
				if (pr.proc) cat("replace cls.", k2, "to", k1,"\n")
				xcl$cluster <- replace(xcl$cluster, (xcl$cluster == k2), k1)
				xcl$size[k1] <- xcl$size[k1] + xcl$size[k2]
				xcl$size[k2] <- 0
				xcl$lnL[k1] <- z$lnLmgd
				xcl$lnL[k2] <- 0
				xcl$detVx[k1] <- z$detVxmgd
				xcl$detVx[k2] <- 0
				xcl$centers[k1,] <- z$ctrmgd
				xcl$centers[k2,] <- 0
				
			}
		}
	}
	list(cluster = xcl$cluster, centers = xcl$centers, size = xcl$size)
}



