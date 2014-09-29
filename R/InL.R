# log-likelihood under the assumption of 
# p-dimensional multivariate normal distribution.
# ignore.covar: ignore the covariance 
# Author: tunenori
###############################################################################

lnL <- function(x, centers, ignore.covar=T){
	x <- as.matrix(x)
	p <- ncol(x)	# p-dimensional multivariate
	n <- nrow(x)	# sample size
	if (missing(centers)) 
		stop("centers must be a number or a matrix")
	if (n <= 2)	# few data
		return(list(lnL=NA, detVx=NA))
	vx <- var(x)	# var-co.var matrix
	# print(x)	# for debug
	if (p == 1){ # x is vector 
		invVx <- 1 / as.vector(vx)
		detVx <- as.vector(vx)
	}else{ 
		if (ignore.covar){
			invVx <- diag(1/diag(vx)) # inv. matrix when assuming diag.  
			detVx <- prod(diag(vx)) # det. when assuming diag. 
		}else{
			invVx <- solve(vx) # inverse matrix of "vx"
			y <- chol(vx) # Cholesky decomposition
			detVx <- prod(diag(y)) # vx = t(y) %*% y, where y is triangular,
			# then, det(vx) = det(t(y)) * det(y)
		}
	}
	t1 <- -p/2 * 1.837877066 # 1.837... = log(2 * 3.1415...)
	t2 <- -log(detVx) / 2
	xmu <- t(apply(x, 1, "-", centers))
	# print(centers)	# for debug
	# print(xmu)	# for debug
	# s <- 0
	# for (i in 1:n)
	#	s <- s + t(xmu[i,]) %*% invVx %*% xmu[i,]
	if (p == 1){
		s <- sum(xmu^2 * invVx)
	}else{
		s <- sum(apply(xmu, 1, txInvVxX, invVx=invVx))
	}
	t3 <- -s / 2
	ll <- (t1 + t2) * n + as.numeric(t3)	# log likelihood
	list(lnL=ll, detVx=detVx)
}

