# A wrapper which combines the clustering algorithm with the continuous distribution calculator	
###############################################################################
tsrFind <- function(TSS_data,binWidth=10,count=1,cluster_Num=1,Draw=FALSE,print=FALSE, countMin=25) {
	TSS_data[,2] -> genes_vector
		length(genes_vector) -> n_genes
		unique(genes_vector) -> unique_genes
		length(unique_genes) -> len_unique
		array(NA,c(len_unique,2)) -> gene_count_array
		match(genes_vector,unique_genes) -> match_array
		print("Calculating the number of TSS per gene in the dataset")
			for (i in 1:len_unique) {
				print(i)
				1:len_unique -> len_array
				len_array[i] -> this_gene
				which(match_array==this_gene) -> gene_index
				length(gene_index) -> gene_count_array[i,2]
				unique_genes[i] -> gene_count_array[i,1]
				}
			as.numeric(gene_count_array[,2]) -> gene_counts
			which(gene_counts<countMin) -> minIndex
			#print(head(minIndex))
			length(minIndex)
			unique_genes[minIndex] -> remove_names
			#print(remove_names)
			match(genes_vector,remove_names) -> remove_index #
			#print(head(remove_index))
			#print(length(remove_index))
			which(is.na(remove_index)) -> keep_index
		TSS_data[keep_index,] -> TSS_data_i #keeping the rows from TSS_data that are more than the specified countMin value
		print("Genes with less than 'countMin' removed from the dataset")

TSR_packagev2 <- function(TSS_data,gene=1,clusterNum=1,binWidth=10,count=1,draw=TRUE) {
	library('moments')
	options(scipen=999) 

	 var_calc_i <- function(x, gene, iterations=1000000) {
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
	cluster_PDF <- function(x,cluster=clusterNum,binwidth=5) {
		if (is(x)!="clustering") {
			stop("Object must be of class 'clustering'")
		}
		total.clusters <- x$clust.string #getting the numbered list of clusters for each TSS tag within the gene set
		unique.clusters <- unique(total.clusters) #listing the total clusters within the gene set ##not sure if this command is necessary
		c_size <- x$cluster.size #provides an array with the list of the number of tags in each cluster
		sort(c_size,decreasing=TRUE) -> cluster_sort #sorts c_size from largest to smallest number of tags in each cluster
		#print(cluster_sort)
		cluster_sort[cluster] -> cluster_pos #selecting the cluster with the largest number of tags
		#print(cluster_pos)
		which(c_size==cluster_pos) -> this_TSS_pos #identifying which of the clusters within the sample is largest (and therefore equal to 'cluster_pos'
		#print(this_TSS_pos) #for debugging
		#print(total.clusters) #for debugging
			if (length(this_TSS_pos) > 1) {
				#print("Top two clusters are of equal size") #for debugging
				length(this_TSS_pos) -> same_len # determines the length of the vector
				#print(same_len)
				x$table -> vc_table #creating a table with cluster numbers and their sizes
				match(this_TSS_pos,vc_table[2,]) -> TSS_table_index
				#print(this_TSS_pos) #for debugging
				#print(vc_table) #for debugging
				#print(TSS_table_index) #for debugging
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
			cluster_array[this_min_index,1] -> this_cluster	#this is the part that I need to focus on
			which(total.clusters==this_cluster) -> this_clust_index
			x$tss[this_clust_index] -> TSSs_pos 
			#print(TSSs_pos) #for debugging
			this_cluster -> this_TSS_pos
 			}
		which(total.clusters == this_TSS_pos) -> TSSs_pos #selecting the positions of the tags that belong to the largest cluster #there seems to be a problem with this line. Began debugging on 2/5/13 at 14:35
		x$tss[TSSs_pos] -> plot_TSSs #selecting the actual tags that belong to the largest cluster from within the var_calc object
		max(plot_TSSs) -> TSS_max #calculating the maximum TSS in the cluster
		min(plot_TSSs) -> TSS_min #caclulating the minimum TSS in the cluster
		abs(TSS_max-TSS_min) -> h_range #calculating the range value between the min and max TSS within the cluster
		#print(h_range)
		h_range/binwidth -> break_n #defining the number of breaks to be in the histogram 
		#print(binwidth)
		if ((break_n > binwidth) == TRUE) { #the regular condition, whereby the number of breaks is larger than the binwidth
			hist(plot_TSSs,breaks=break_n,plot=FALSE) -> output_hist 
		}
		if ((break_n <= binwidth) == TRUE) { ##this is a problem area. Please look into this further #2.5.13
			binwidth -> break_n #not sure why we need this #revisit
			hist(plot_TSSs,breaks=break_n,plot=FALSE) -> output_hist
		}
		#hist(plot_TSSs,xlim=c(TSS_min,TSS_max),breaks=break_n,plot=FALSE) -> output_hist
		#lines(density(plot_TSSs,kernel="gaussian"), col='red', lwd=3)
		dimnames(x$table)[[1]][1] -> gene_name
		length(plot_TSSs) -> n_TSS
		array(NA,c(n_TSS,1)) -> TSS_array
		plot_TSSs -> TSS_array[,1]
		gene_name -> colnames(TSS_array) 
		as.data.frame(TSS_array) -> TSS_array
		list(table=plot_TSSs,hist=output_hist,array=TSS_array) -> TSS_List
		return(TSS_List)
	}
	hist_peak_find <- function(x,count=1) {
	hist_obj <- x$hist # the 'hist' object that we will be retrieving data from
	#print(hist_obj) #for debugging
	hist_table <- x$table
	hist_counts <- hist_obj$counts # the 'counts' vector that our algorithm is focused on analyzing
	hist_mids <- hist_obj$mids # the 'mids' vector that our algorithm is focused on analyzing
	hist_cont <- which(hist_counts>=count) # an index of the positions within 'hist_counts' for which there are 'count' or more TSS tags
	ext_counts <- hist_counts[hist_cont]
	#print(ext_counts) #for debugging
	mid_extr <- hist_mids[hist_cont] 
	cont_len <- length(hist_cont) #the length of the vector 'hist_cont', which we'll use in the subsequent 'for' loop
	#print(cont_len) #for debugging
	#print(ext_counts) #for debugging
	cont_array <- array(NA,c((cont_len),3))
	cont_array[,1] <- 1:(cont_len)
	#print(hist_cont) #for debugging
	if (length(hist_mids)>1) {
		if (cont_len==length(hist_counts)) #if all of the values have at least 'count' in them
		{
			cont_len -> cont_array[,2]
		}
		else {
			if (hist_cont[2]-hist_cont[1]>1) #starting the loop for i=1
			{ 
				1 -> end_pos 
				end_pos -> cont_array[1,2]
			}
			for (i in 2:cont_len) { #continuing the loop 
				i -> start_pos
				#if (i<cont_len){
				i -> end_pos 
				#}
				if ((hist_cont[i]-hist_cont[i-1])>1) { #if i and i-1 are not adjacent
					start_pos -> end_pos
					end_pos -> cont_array[i,2]
					next
				}
				if ((hist_cont[i]-hist_cont[i-1])==1) { #if i and i-1 ARE adjacent (the contiguous strings we are looking for)
					i -> this_start
					cont_len-this_start -> remain_len #defines how much longer the loop will continue for
					#	print(remain_len) #for debugging
					for (j in 1:remain_len) { #loop within
						if ((i+j)>cont_len) { break }
						if ( hist_cont[i+j]-hist_cont[i+(j-1)]==1 ){ #flow control issues- gives NA/NA in last two columns of last row when evaluated
							i+j -> end_pos #added for testing
							end_pos -> cont_array[this_start:i,2] #added for testing
							next
						}
						if ((hist_cont[i+j]-hist_cont[i+j-1])>1) { 
							i+j -> end_pos
							end_pos -> cont_array[this_start:i,2] #added 'this_start' to create a top boundary
							break
						}
					}
					(i + j) -> i
				}
				NULL -> this_start
				NULL -> start_pos
				NULL -> end_pos
			}
			cont_array[1,2] <- cont_array[1,1]	
			cont_array[cont_len,2] <- cont_array[cont_len,1]
		}
		if (cont_len==1){
			max(ext_counts) -> max_ext
			which(ext_counts==max_ext) -> max_pos
			ext_counts[max_pos] -> cont_values
			mid_extr[max_pos] -> mids_values
			#print(hist_cont)
			#print(max_ext)
			#print(ext_counts)	
			#print(max_pos)
			#print(mids_values)
			#print(cont_values)
		}
		if (cont_len>1){ ##this is where I neeed to continue debuggging 5/1/12
			if (hist_cont[cont_len]-hist_cont[cont_len-1]>1) { #a loop to remove the last row of cont_array if the two values are not contiguous
				cont_array[-cont_len,] -> cont_array	
			}
			if (is.array(cont_array)){
				#print(cont_array) for debugging
				for (k in 1:length(cont_array[,1])){
					cont_array[k,2]-cont_array[k,1] -> cont_array[k,3]
				}
			}
			else {
				cont_array[1]-cont_array[2] -> cont_array[3]
			}
			if (is.array(cont_array)) {
				max_value <- max(cont_array[,3])
				#print(max_value)
				max_pos <- which(cont_array[,3]==max_value)
				if (length(max_pos)==1&&cont_array[max_pos,3]>2){ #modified the algorithm to include very peaked TSRs #done on 5/1/12
					#print(max_pos)
					#print(cont_array)
					#print(cont_len)
					#print(hist_cont)
					#print(ext_counts)
					#print(mid_extr)
					#print(cont_array)
					max_coord <- cont_array[max_pos,1:2]
					cont_values <- ext_counts[max_coord[1]:max_coord[2]] #
					#mids_values <- hist_mids[max_coord[1]:max_coord[2]]
					mids_values <- mid_extr[max_coord[1]:max_coord[2]]
					if (length(hist_cont)>1){
						if (hist_cont[cont_len]-hist_cont[cont_len-1]>1) { #another loop to remove the last row of cont_array if the two values are not contiguous
							max_coord[2] <- (cont_array[max_pos,2])-1
							cont_values <- ext_counts[max_coord[1]:max_coord[2]] #
							mids_values <- mid_extr[max_coord[1]:max_coord[2]]
						}
					}
				}
				else {
					max(ext_counts) -> max_ext
					which(ext_counts==max_ext) -> max_pos
					ext_counts[max_pos] -> cont_values
					mid_extr[max_pos] -> mids_values
					#print("Row 118")
					#print(max_pos)
				}
			}
			else {
				max(ext_counts) -> max_ext
				which(ext_counts==max_ext) -> max_pos
				ext_counts[max_pos] -> cont_values
				mid_extr[max_pos] -> mids_values
				#print("Row 127")
				#print(max_ext)
				#print(cont_array)
			}
		}
	}
	if (length(hist_mids)==1){
		cont_values <- length(hist_table) 
		mids_values <- mean(hist_table)
	}
	#print(mid_extr) #for debugging
	#print(max_coord) #for debugging
	#print(cont_values) #for debugging
	#print(length(cont_values)) #for debugging
	#print(length(max_coord[1]:max_coord[2])) #for debugging
	hist_array <- array(NA,c(length(cont_values),2))
	hist_array[,1] <- cont_values
	hist_array[,2] <- mids_values
	colnames(hist_array) <- c("bin_counts","bin_midpoints") #naming the columns for the array
	hist_return <- list(array=hist_array,vector=hist_table)
	return(hist_return)
}
	PDF_fit <- function(x,normal=TRUE) {
		h <- hist(x) #changed for this function
		xhist <- c(min(h$breaks),h$breaks)
		yhist <- c(0,h$density,0)
		xfit <- seq(min(x),max(x),length=40) # changed from original to fit this package
		if (normal==TRUE){
			yfit <- dnorm(xfit,mean=mean(x),sd=sd(x)) #also changed from original to fit this package
		}
		if (normal==FALSE) {
			yfit <- dlogis(xfit,location=mean(x),scale=sd(x)) #also changed from original to fit this package
		}
		plot(xhist,yhist,type="s",ylim=c(0,max(yhist,yfit)),main="Normal pdf and histogram")
		lines(xfit,yfit,col='red2')
	}
	TSS_clustered <- var_calc_i(TSS_data_i,gene) #creating a clustering object using var_calc_i
	PDF_object <- cluster_PDF(TSS_clustered,binwidth=binWidth) #isolating the major cluster of TSSs according to PDF_object, with a bidwidth of the selected size
	#print(str(PDF_object))
	hist_peak <- hist_peak_find(PDF_object,count) #detect the longest continuously occupied stretch of bins in the cluster
	#print(hist_peak)
	segment_vector2 <- hist_peak$vector #inserting the calculated contiguous string of values
	vec_range <- range(hist_peak$array[,2])
	if (length(vec_range>1)){
	vec_min <- vec_range[1]-(binWidth/2)
	vec_max <- vec_range[2]+(binWidth/2)
}
	else {
		vec_min <- hist_peak$array[1,2]-(binWidth/2)
		vec_min <- hist_peak$array[1,2]+(binWidth/2)
	}
	segment_index <- which(segment_vector2 >= vec_min & segment_vector2 <= vec_max )
	segment_vector3 <- segment_vector2[segment_index]
	#print(segment_vector2)
	#print(segment_vector3)
	#print(hist_peak$array)
	#print(hist_vector)
	#hist_Array <- hist_peak$array
	#hist_array_len <- length(hist_Array[,1])
	#segment_min <- hist_Array[1,2]
	#segment_max <- hist_Array[hist_array_len,2]
	#print(segment_min)
	#print(segment_max)
	#segment_vector <- hist_vector[(hist_vector > segment_min) & (hist_vector < segment_max)] 
	#hist_seg1 <- hist_vector[segment_1]
	#segment_2 <- which(hist_seg1<=segment_max)
	#segment_vector <- hist_vector[segment_index]
	TSR_shape <- kurtosis(segment_vector3)
	if(is.na(TSR_shape)) {
		c('Completely_Peaked') -> TSR_shape
	}
	TSR_breadth <- max(segment_vector3)-min(segment_vector3)+1
	TSR_mid <- mean(segment_vector3)
	TSR_mid <- round(TSR_mid,0)
	TSR_range <- range(segment_vector3)
	nTags <- length(segment_vector3)
	total_tags <- length(TSS_clustered$tss)
	peak_ratio <- nTags/total_tags
	TSR_output <- list(shape=TSR_shape,breadth=TSR_breadth,range=TSR_range,mid=TSR_mid,tss_vector=segment_vector3,tag_count=nTags,total_count=total_tags,Ratio=peak_ratio)
	if (draw==TRUE) {
		PDF_fit(segment_vector2,normal=TRUE)
	}
	return(TSR_output)
}
	TSS_data_i[,2] -> genes_vector
	unique(genes_vector) -> genesUnique
	length(genesUnique) -> ngenes
	TSR_array <- array(NA,c(ngenes,8))
	colnames(TSR_array) <- c('Shape','TSR_Breadth','Peak_Count','Total_Count','Tag_Ratio','Midpoint','Start','End')
	rownames(TSR_array) <- genesUnique
	for (t in 1:ngenes) {
		genesUnique[t] -> this_Gene
		print(this_Gene)
		TSR_packagev2(TSS_data,this_Gene,cluster_Num,binWidth,count,Draw) -> TSR_out
		TSR_out$shape -> TSR_array[t,1]
		TSR_out$breadth -> TSR_array[t,2]
		TSR_out$tag_count -> TSR_array[t,3]
		TSR_out$total_count -> TSR_array[t,4]
		TSR_out$Ratio -> TSR_array[t,5]
		TSR_out$mid -> TSR_array[t,6]
		TSR_out$range -> TSR_array[t,7:8]
	}
	if (print) {
		write.table(TSR_array,file="TSR_output.txt",sep="\t",row.names=TRUE,col.names=TRUE)
	}
	return(TSR_array)
}
