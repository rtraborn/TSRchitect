# Xmeans Script from U of C
# non-hierarchal clustering
# $Id: xmeans.prog,v 1.21 2010/02/03 01:44:06 tunenori Exp tunenori $
# n: sample size 
# p: dimension of data
# it: iteration
# Author: tunenori
###############################################################################

xmeansMain <- function(n = 100, p = 2, it = 1000){
	#library(mva)
	library(stats)
	
	ik <- 2 # initial cluster number; default = 2
	iter.max <- 20 # maximum of iterations applied for "kmeans"
	
	#x<-rbind(matrix(rnorm(n,mean=0,sd=0.2),ncol=p),
	#matrix(rnorm(n,mean=1,sd=0.2),ncol=p),
	#matrix(rnorm(n,mean=2,sd=0.2),ncol=p),
	#matrix(rnorm(n,mean=3,sd=0.2),ncol=p),
	#matrix(rnorm(n,mean=-1,sd=0.2),ncol=p))
	
	np <- n / p
	x<-rbind(cbind(rnorm(2*np,mean=0,sd=0.2),rnorm(np,mean=0,sd=0.2)),
			cbind(rnorm(np,mean=-2,sd=0.3),rnorm(np,mean=0,sd=0.3)),
			cbind(rnorm(np,mean=2,sd=0.3),rnorm(np,mean=0,sd=0.3)),
			cbind(rnorm(np,mean=0,sd=0.4),rnorm(np,mean=2,sd=0.4)),
			cbind(rnorm(np,mean=0,sd=0.4),rnorm(np,mean=-2,sd=0.4)))
	r <- 0.0
	r1 <- sqrt(1-r^2)
	#postscript("xplotr05landscape.ps", horizontal=F, width=5, height=5.5)
	#plot(x)
	#graphics.off()
	
	n <- NULL
	ne <- NULL # effecive cluster number list
	size <- list()
	for (i in 1:it){
		
		u0 <- rnorm(2*np)
		u1 <- r * u0 + r1 * rnorm(2*np)
		c1 <- cbind(0 + 0.2 * u0, 0 + 0.2 * u1) 
		# c1 <- cbind(0 + 0.2 * u0, 0 + 0.2 * u1) 
		u0 <- rnorm(np)
		u1 <- r * u0 + r1 * rnorm(np)
		c2 <- cbind(-2 + 0.3 * u0, 0 + 0.3 * u1) 
		# c2 <- cbind(-1 + 0.2 * u0, -1 + 0.2 * u1) 
		u0 <- rnorm(np)
		u1 <- r * u0 + r1 * rnorm(np)
		c3 <- cbind(2 + 0.3 * u0, 0 + 0.3 * u1) 
		# c3 <- cbind(1 + 0.2 * u0, 1 + 0.2 * u1) 
		u0 <- rnorm(np)
		u1 <- r * u0 + r1 * rnorm(np)
		c4 <- cbind(0 + 0.4 * u0, 2 + 0.4 * u1) 
		# c4 <- cbind(2 + 0.2 * u0, 2 + 0.2 * u1) 
		u0 <- rnorm(np)
		u1 <- r * u0 + r1 * rnorm(np)
		c5 <- cbind(0 + 0.4 * u0, -2 + 0.4 * u1) 
		# c5 <- cbind(3 + 0.2 * u0, 3 + 0.2 * u1) 
		x<-rbind(c1, c2, c3, c4, c5)
		
		xcl <- xmeans(x, ik, iter.max, merge.cls=T, pr.proc = T)
		plot(x,type="n")
		text(x, labels=xcl$cluster)
		cat(paste("i =", i, "\n"))
		cat("xcl$size =")
		print(xcl$size)
		n <- c(n, length(xcl$size))
		ne <- c(ne, length(xcl$size[xcl$size != 0]))
		if (length(xcl$size[xcl$size != 0]) == 1) break
		size <- c(size, list(xcl$size))
	}
	
	cat("table(ne):\n")
	table(ne)
}
