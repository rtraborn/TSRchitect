tsrLook  <- function(TSS_data,gene=1,clusterNum=1,binWidth=10,count=1,draw=TRUE) {
	library('moments')
	options(scipen=999) 
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
	TSS_clustered <- varCalc(TSS_data_i,gene) #creating a clustering object using varCalc
	PDF_object <- cluster_PDF(TSS_clustered,binwidth=binWidth) #isolating the major cluster of TSSs according to PDF_object, with a bidwidth of the selected size
	hist_peak <- histFind(PDF_object,count) #detect the longest continuously occupied stretch of bins in the cluster
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