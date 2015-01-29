varCalc <- function(x, gene, iterations=1000000) {
 		gene -> gene_name
 		which(x[,2]==gene) -> this_gene_index
 		x[this_gene_index,3:5] -> tss_slice
 		tss_slice[1,1] -> this_strand
 		if (this_strand=="+") {
 		tss_slice[,2] -> tss_vector 
 		}
 		else {tss_slice[,3] -> tss_vector}
 		unique(tss_vector) -> unique_tss
 		length(unique_tss) -> this_unique
 			if (this_unique < 2) {
 				array(NA,c(7,1)) -> table_array
 				unique_tss -> this_centers
 				1 -> this_clusters
 				length(tss_vector) -> this_sizes
 				tss_vector -> this_tss
 				c(NA,1,this_centers,0,this_sizes,100,length(tss_vector)) -> table_array[,1]
 				}
 			else {
 		xmeans(tss_vector,ik = 2, iter.max = iterations, pr.proc = F, ignore.covar=T, merge.cls=F) -> xmeans_output
 		xmeans_output$centers -> this_centers
 		xmeans_output$cluster -> this_clusters
 		xmeans_output$centers -> this_centers
 		xmeans_output$size -> this_sizes
 		tss_vector -> this_tss
		
		nrow(this_centers) -> nclusters 	#creating the clusters_array
		array(NA,c(3,nclusters)) -> clusters_array 
		1:nclusters -> clusters_array[1,]
		this_centers[,1] -> clusters_array[2,]
		this_sizes -> clusters_array[3,]
		c("Cluster Number","Center","Size") -> rownames(clusters_array)
		
		array(NA,c(7,nclusters)) -> table_array
		1:nclusters -> table_array[2,]
		this_centers -> table_array[3,]
		this_sizes -> table_array[5,]
		sum(this_sizes) -> table_array[7,]
		table_array[7,1] -> total_tss
			for (i in 1:nclusters) {
			which(this_clusters==i) -> cluster_index
			this_tss[cluster_index] -> cluster_tss
			var(cluster_tss) -> table_array[4,i] 
			(this_sizes[i]/total_tss) -> table_array[6,i]
			}
			}
		rownames(table_array) <- c(gene_name,"Cluster Number","Cluster Centers","Cluster Variance","Cluster Size","Percentage of Total","Total Number of Clones") 
		var_out <- list(table=table_array,centers=this_centers,clusters=this_clusters,clust.string=this_clusters,tss=this_tss,variance=table_array[4,],cluster.size=this_sizes)
		class(var_out) = "clustering"
		return(var_out)
		}