cluster_PDF <- function(x,cluster=clusterNum,binwidth=5) {
		if (is(x)!="clustering") {
			stop("Object must be of class 'clustering'")
		}
		total.clusters <- x$clust.string #getting the numbered list of clusters for each TSS tag within the gene set
		unique.clusters <- unique(total.clusters) #listing the total clusters within the gene set ##not sure if this command is necessary
		c_size <- x$cluster.size #provides an array with the list of the number of tags in each cluster
		sort(c_size,decreasing=TRUE) -> cluster_sort #sorts c_size from largest to smallest number of tags in each cluster
		cluster_sort[cluster] -> cluster_pos #selecting the cluster with the largest number of tags
		which(c_size==cluster_pos) -> this_TSS_pos #identifying which of the clusters within the sample is largest (and therefore equal to 'cluster_pos'
		#print(this_TSS_pos) #for debugging
		#print(total.clusters) #for debugging
			if (length(this_TSS_pos) > 1) {
				#print("Top two clusters are of equal size") #for debugging
				length(this_TSS_pos) -> same_len # determines the length of the vector
				x$table -> vc_table #creating a table with cluster numbers and their sizes
				match(this_TSS_pos,vc_table[2,]) -> TSS_table_index
				array(NA,c(same_len,2)) -> cluster_array
					for (i in 1:same_len) { #loop for as many clusters as are in this vector. (The most likely case is 2, but length > 2 instances could theoretically exist also.)
						which(total.clusters==TSS_table_index[i]) -> cluster_index
						#print(cluster_index) #for debugging
						x$tss[cluster_index] -> this_cluster_array
						#print(this_cluster_array) #for debugging
						range(this_cluster_array) -> this_range1
						abs(this_range1[2]-this_range1[1]) -> this_range2
						TSS_table_index[i]-> cluster_array[i,1] #inserts the cluster number in the first column of the array
						this_range2 -> cluster_array[i,2]
						#print(cluster_array)
					}
			min(cluster_array[,2]) -> this_min
			which(cluster_array[,2]==this_min) -> this_min_index
			cluster_array[this_min_index,1] -> this_cluster	
			which(total.clusters==this_cluster) -> this_clust_index
			x$tss[this_clust_index] -> TSSs_pos 
			#print(TSSs_pos) #for debugging
			this_cluster -> this_TSS_pos
 			}
		which(total.clusters == this_TSS_pos) -> TSSs_pos #selecting the positions of the tags that belong to the largest cluster
		x$tss[TSSs_pos] -> plot_TSSs #selecting the actual tags that belong to the largest cluster from within the var_calc object
		max(plot_TSSs) -> TSS_max #calculating the maximum TSS in the cluster
		min(plot_TSSs) -> TSS_min #caclulating the minimum TSS in the cluster
		abs(TSS_max-TSS_min) -> h_range #calculating the range value between the min and max TSS within the cluster
		h_range/binwidth -> break_n #defining the number of breaks to be in the histogram 
		if ((break_n > binwidth) == TRUE) { #the regular condition, whereby the number of breaks is larger than the binwidth
			hist(plot_TSSs,breaks=break_n,plot=FALSE) -> output_hist 
		}
		if ((break_n <= binwidth) == TRUE) { ##this is a problem area. Please look into this further #2.5.13
			binwidth -> break_n #not sure why we need this #revisit
			hist(plot_TSSs,breaks=break_n,plot=FALSE) -> output_hist
		}
		dimnames(x$table)[[1]][1] -> gene_name
		length(plot_TSSs) -> n_TSS
		array(NA,c(n_TSS,1)) -> TSS_array
		plot_TSSs -> TSS_array[,1]
		gene_name -> colnames(TSS_array) 
		as.data.frame(TSS_array) -> TSS_array
		list(table=plot_TSSs,hist=output_hist,array=TSS_array) -> TSS_List
		return(TSS_List)
	}