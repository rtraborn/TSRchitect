counts_filter <- function(TSS_data, countMin=25) {
	TSS_data[,2] -> genes_vector
		length(genes_vector) -> n_genes
		unique(genes_vector) -> unique_genes
		length(unique_genes) -> len_unique
		array(NA,c(len_unique,2)) -> gene_count_array
		match(genes_vector,unique_genes) -> match_array
		print("Calculating the number of TSS per gene in the dataset")
			for (i in 1:len_unique) {
				#print(i)
				1:len_unique -> len_array
				len_array[i] -> this_gene
				which(match_array==this_gene) -> gene_index
				length(gene_index) -> gene_count_array[i,2]
				print(length(gene_index))
				unique_genes[i] -> gene_count_array[i,1]
				}
			as.numeric(gene_count_array[,2]) -> gene_counts
			which(gene_counts<countMin) -> minIndex
			#print(head(minIndex))
			length(minIndex)
			unique_genes[minIndex] -> remove_names
			print(remove_names)
			match(genes_vector,remove_names) -> remove_index #need to fix this line- need a set operation than permits duplicate values and overapping
			#print(head(remove_index))
			print(length(remove_index))
			which(is.na(remove_index)) -> keep_index
		TSS_data[keep_index,] -> TSS_data_i #keeping the rows from TSS_data that are more than the specified countMin value
		print("Genes with less than 'countMin' removed from the dataset")
		return(TSS_data_i)
		}