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
			length(minIndex)
			unique_genes[minIndex] -> remove_names
			match(genes_vector,remove_names) -> remove_index #
			which(is.na(remove_index)) -> keep_index
		TSS_data[keep_index,] -> TSS_data_i #keeping the rows from TSS_data that are more than the specified countMin value
	TSS_data_i[,2] -> genes_vector
	unique(genes_vector) -> genesUnique
	length(genesUnique) -> ngenes
	TSR_array <- array(NA,c(ngenes,8))
	colnames(TSR_array) <- c('Shape','TSR_Breadth','Peak_Count','Total_Count','Tag_Ratio','Midpoint','Start','End')
	rownames(TSR_array) <- genesUnique
	for (t in 1:ngenes) {
		genesUnique[t] -> this_Gene
		print(this_Gene)
		tsrLook(TSS_data,this_Gene,cluster_Num,binWidth,count,Draw) -> TSR_out
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
