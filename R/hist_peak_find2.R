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