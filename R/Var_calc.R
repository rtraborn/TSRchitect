# works well as of 2.15.2011 at 4:15PM #slightly adjusted on 2.18.2011 at 7:23PM #$cluster output modified slightly on 3.2.2011 at 4:08PM
# Calculates the variance within all clusters using the outcome from 'x_means_mod_2'
# Author: R. Taylor Raborn
###############################################################################

var_calc <- function(x, gene=1,iterations=1000000) {
	xcenters<- function(x, gene=1, center=2, iterations = 100000) {
		onecenter <- function(x, gene=1) { #function that specifies which column (gene) to select
			z <- array(c(na.omit(x[,gene])), dim=c(length(na.omit(x[,gene])),1))
			x1 <- as.vector(z[1])
			x1 <- array(x1,dim=c(1,1))
			xx <- array(NA,dim=c(7,1)) #creates an array of the appropriate size- in this case 1.
			xx[2,] <- c("1")
			xx[3,] <- x1
			xx[4,] <- c("0") #filling in the new array with the center = 1 data
			xx[5,] <- length(z)
		    xx[6,] <- c("100") 
			xx[7,] <- length(z) #entering the amount of clones for this gene
			return(xx)
		}
		z <- array(c(na.omit(x[,gene])), dim=c(length(na.omit(x[,gene])),1))
		u <- length(unique(z))
		if (u > 1) {
		list.i <- xmeans_mod(z, ik = center, iter.max = iterations, pr.proc = F, ignore.covar=T, merge.cls=F)
		x1 <- round(list.i$centers,2) #using the xmeans modified function
		x1 <- array(x1,dim=c(1,length(x1))) #turning the above into an array of the apporiate dimensions
		x2 <- list.i$size #using the kmeans modified function
		x2 <- array(x2,dim=c(1,length(x2)))
		x3 <- list.i$cluster # a list of positions of the clusters in the array
		x4 <- as.vector(z) # a list of the data points that correspond to the clusters in x3
		x5 <- list.i$size
		if (length(x2) > length(x1)) { x1 <- c(x1,rep(NA,c(length(x2)-length(x1)))) }#making sure the dimensions match for both arrays
		if (length(x1) > length(x2)) { x2 <- c(x2,rep(NA,c(length(x1)-length(x2)))) }
		xx <- array(NA,dim=c(7,length(x1))) #creates an array of the appropriate size
		xx[3,] <- x1 #filling in the new array with the xmeans data, part 1
		xx[5,] <- x2 #filling in the new array with the xmeans data, part 2
		xx[6,] <- NA #making a blank row, to be filled with the '% of total' later
		xx[7,] <- length(z) #entering the amount of clones for this gene
		xcen <- list(table=xx,clusters=x3, centers=x1, tss=x4, size=x5)
	}
	else {
		xcen <- onecenter(x,gene)
		xcen <- list(table=xcen,clusters=xcen[2,],centers=xcen[3,],tss=as.vector(z),size=as.vector(xcen[5,]))
	}
  return(xcen)  	
}
	tss.table <- xcenters(x,gene,iterations=iterations) #generating the table to be called later in this program
	cluster.list <- tss.table$clusters #generates an array with the positions of the clusters
	tss.array <- tss.table$tss #an array with the TSSs for that gene
	if (length(tss.table$size) == 1) {
		var.array <- tss.table$table
		var.array[2,1] <- 1
		var.array[6,1] <- 100
		var.i <- 0
	}
	else {
	tss.total <- sum(tss.table$size)
	max.clust <- max(tss.table$size)
	clust.num <- length(tss.table$size) #finding the number of distinct clusters in the array
	cluster.array <- array(NA, c(max.clust,clust.num))	
	var.array <- array(NA, c(7,clust.num))
		for (j in 1:clust.num) {
			same.clust <- which(cluster.list == j) 
			in.clust <- length(same.clust) #how many are in this cluster?
			for (i in 1:in.clust) {
				cluster.array[i,j] <- tss.array[same.clust[i]]
			}
		}
	for (i in 1:clust.num) {
		tss.list <- na.omit(cluster.array[,i]) 
		tss.list <- as.vector(tss.list)
		var <- var(tss.list)
		var <- round(var,2)
		var.array[4,i] <- var
	}
		var.array[1,] <- rep(NA,clust.num)
		var.array[2,] <- 1:clust.num
		var.array[3,] <- tss.table$centers
		var.array[5,] <- tss.table$size
		var.array[7,] <- tss.total
		for (i in 1:clust.num) {
			c.size <- as.numeric(var.array[5,i])
			t.size <- as.numeric(var.array[7,i])
			clust.per <- (c.size/t.size)*100
			var.array[6,i] <- clust.per
		}
		var.i <- var.array[4,]
	}
		gene_names <- colnames(x)
		gene_name <- gene_names[gene]
		clust.array <- unique(cluster.list)
		center.clusters <- array(NA,c(3,length(clust.array)))
		center.clusters[1,] <- var.array[2,]
		center.clusters[2,] <- var.array[3,]
		center.clusters[3,] <- var.array[5,]
		rownames(center.clusters) <- c("Cluster Number","Center","Size")
	rownames(var.array) <- c(gene_name,"Cluster Number","Cluster Centers","Cluster Variance","Cluster Size","Percentage of Total","Total Number of Clones")
	var_out <- list(table=var.array,centers=var.array[3,],clusters=center.clusters,clust.string=cluster.list, tss=tss.array,variance=var.i,cluster.size=as.numeric(tss.table$size))
	class(var_out) = "clustering"
	return(var_out)
	}
				
#this currently generates the following:

#> var_calc(miura_rich,2,iterations=1000000)$table
#                            [,1]       [,2]
#YFL014W                       NA         NA
#Cluster Number           2.00000   1.000000
#Cluster Centers         59.27000  27.600000
#Cluster Variance         7.00000 109.000000
#Cluster Size           391.00000   5.000000
#Percentage of Total     98.73737   1.262626
#Total Number of Clones 396.00000 396.000000

		
		
	

