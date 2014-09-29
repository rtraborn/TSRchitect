\name{TSR_Completev2}
\alias{TSR_Completev2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Identifies the putative promoter (TSR) from each expressed gene
  within a large-scale 5' end dataset
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
TSR_Completev2(TSS_data, binWidth = 10, count = 1, cluster_Num = 1, Draw = FALSE, print = FALSE, countMin = 25)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{TSS_data}{
%%     ~~Describe \code{TSS_data} here~~
}
  \item{binWidth}{
%%     ~~Describe \code{binWidth} here~~
}
  \item{count}{
%%     ~~Describe \code{count} here~~
}
  \item{cluster_Num}{
%%     ~~Describe \code{cluster_Num} here~~
}
  \item{Draw}{
%%     ~~Describe \code{Draw} here~~
}
  \item{print}{
%%     ~~Describe \code{print} here~~
}
  \item{countMin}{
%%     ~~Describe \code{countMin} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (TSS_data, binWidth = 10, count = 1, cluster_Num = 1, 
    Draw = FALSE, print = FALSE, countMin = 25) 
{
    genes_vector <- TSS_data[, 2]
    n_genes <- length(genes_vector)
    unique_genes <- unique(genes_vector)
    len_unique <- length(unique_genes)
    gene_count_array <- array(NA, c(len_unique, 2))
    match_array <- match(genes_vector, unique_genes)
    print("Calculating the number of TSS per gene in the dataset")
    for (i in 1:len_unique) {
        print(i)
        len_array <- 1:len_unique
        this_gene <- len_array[i]
        gene_index <- which(match_array == this_gene)
        gene_count_array[i, 2] <- length(gene_index)
        gene_count_array[i, 1] <- unique_genes[i]
    }
    gene_counts <- as.numeric(gene_count_array[, 2])
    minIndex <- which(gene_counts < countMin)
    length(minIndex)
    remove_names <- unique_genes[minIndex]
    remove_index <- match(genes_vector, remove_names)
    keep_index <- which(is.na(remove_index))
    TSS_data_i <- TSS_data[keep_index, ]
    print("Genes with less than 'countMin' removed from the dataset")
    TSR_packagev2 <- function(TSS_data, gene = 1, clusterNum = 1, 
        binWidth = 10, count = 1, draw = TRUE) {
        library("moments")
        options(scipen = 999)
        var_calc_i <- function(x, gene, iterations = 1e+06) {
            gene_name <- gene
            this_gene_index <- which(x[, 2] == gene)
            tss_slice <- x[this_gene_index, 3:5]
            this_strand <- tss_slice[1, 1]
            if (this_strand == "+") {
                tss_vector <- tss_slice[, 2]
            }
            else {
                tss_vector <- tss_slice[, 3]
            }
            unique_tss <- unique(tss_vector)
            this_unique <- length(unique_tss)
            if (this_unique < 2) {
                table_array <- array(NA, c(7, 1))
                this_centers <- unique_tss
                this_clusters <- 1
                this_sizes <- length(tss_vector)
                this_tss <- tss_vector
                table_array[, 1] <- c(NA, 1, this_centers, 0, 
                  this_sizes, 100, length(tss_vector))
            }
            else {
                xmeans_output <- xmeans_mod(tss_vector, ik = 2, 
                  iter.max = iterations, pr.proc = F, ignore.covar = T, 
                  merge.cls = F)
                this_centers <- xmeans_output$centers
                this_clusters <- xmeans_output$cluster
                this_centers <- xmeans_output$centers
                this_sizes <- xmeans_output$size
                this_tss <- tss_vector
                nclusters <- nrow(this_centers)
                clusters_array <- array(NA, c(3, nclusters))
                clusters_array[1, ] <- 1:nclusters
                clusters_array[2, ] <- this_centers[, 1]
                clusters_array[3, ] <- this_sizes
                rownames(clusters_array) <- c("Cluster Number", 
                  "Center", "Size")
                table_array <- array(NA, c(7, nclusters))
                table_array[2, ] <- 1:nclusters
                table_array[3, ] <- this_centers
                table_array[5, ] <- this_sizes
                table_array[7, ] <- sum(this_sizes)
                total_tss <- table_array[7, 1]
                for (i in 1:nclusters) {
                  cluster_index <- which(this_clusters == i)
                  cluster_tss <- this_tss[cluster_index]
                  table_array[4, i] <- var(cluster_tss)
                  table_array[6, i] <- (this_sizes[i]/total_tss)
                }
            }
            rownames(table_array) <- c(gene_name, "Cluster Number", 
                "Cluster Centers", "Cluster Variance", "Cluster Size", 
                "Percentage of Total", "Total Number of Clones")
            var_out <- list(table = table_array, centers = this_centers, 
                clusters = this_clusters, clust.string = this_clusters, 
                tss = this_tss, variance = table_array[4, ], 
                cluster.size = this_sizes)
            class(var_out) = "clustering"
            return(var_out)
        }
        cluster_PDF <- function(x, cluster = clusterNum, binwidth = 5) {
            if (is(x) != "clustering") {
                stop("Object must be of class 'clustering'")
            }
            total.clusters <- x$clust.string
            unique.clusters <- unique(total.clusters)
            c_size <- x$cluster.size
            cluster_sort <- sort(c_size, decreasing = TRUE)
            cluster_pos <- cluster_sort[cluster]
            this_TSS_pos <- which(c_size == cluster_pos)
            if (length(this_TSS_pos) > 1) {
                same_len <- length(this_TSS_pos)
                vc_table <- x$table
                TSS_table_index <- match(this_TSS_pos, vc_table[2, 
                  ])
                cluster_array <- array(NA, c(same_len, 2))
                for (i in 1:same_len) {
                  cluster_index <- which(total.clusters == TSS_table_index[i])
                  this_cluster_array <- x$tss[cluster_index]
                  this_range1 <- range(this_cluster_array)
                  this_range2 <- abs(this_range1[2] - this_range1[1])
                  cluster_array[i, 1] <- TSS_table_index[i]
                  cluster_array[i, 2] <- this_range2
                }
                this_min <- min(cluster_array[, 2])
                this_min_index <- which(cluster_array[, 2] == 
                  this_min)
                this_cluster <- cluster_array[this_min_index, 
                  1]
                this_clust_index <- which(total.clusters == this_cluster)
                TSSs_pos <- x$tss[this_clust_index]
                this_TSS_pos <- this_cluster
            }
            TSSs_pos <- which(total.clusters == this_TSS_pos)
            plot_TSSs <- x$tss[TSSs_pos]
            TSS_max <- max(plot_TSSs)
            TSS_min <- min(plot_TSSs)
            h_range <- abs(TSS_max - TSS_min)
            break_n <- h_range/binwidth
            if ((break_n > binwidth) == TRUE) {
                output_hist <- hist(plot_TSSs, breaks = break_n, 
                  plot = FALSE)
            }
            if ((break_n <= binwidth) == TRUE) {
                break_n <- binwidth
                output_hist <- hist(plot_TSSs, breaks = break_n, 
                  plot = FALSE)
            }
            gene_name <- dimnames(x$table)[[1]][1]
            n_TSS <- length(plot_TSSs)
            TSS_array <- array(NA, c(n_TSS, 1))
            TSS_array[, 1] <- plot_TSSs
            colnames(TSS_array) <- gene_name
            TSS_array <- as.data.frame(TSS_array)
            TSS_List <- list(table = plot_TSSs, hist = output_hist, 
                array = TSS_array)
            return(TSS_List)
        }
        hist_peak_find <- function(x, count = 1) {
            hist_obj <- x$hist
            hist_table <- x$table
            hist_counts <- hist_obj$counts
            hist_mids <- hist_obj$mids
            hist_cont <- which(hist_counts >= count)
            ext_counts <- hist_counts[hist_cont]
            mid_extr <- hist_mids[hist_cont]
            cont_len <- length(hist_cont)
            cont_array <- array(NA, c((cont_len), 3))
            cont_array[, 1] <- 1:(cont_len)
            if (length(hist_mids) > 1) {
                if (cont_len == length(hist_counts)) {
                  cont_array[, 2] <- cont_len
                }
                else {
                  if (hist_cont[2] - hist_cont[1] > 1) {
                    end_pos <- 1
                    cont_array[1, 2] <- end_pos
                  }
                  for (i in 2:cont_len) {
                    start_pos <- i
                    end_pos <- i
                    if ((hist_cont[i] - hist_cont[i - 1]) > 1) {
                      end_pos <- start_pos
                      cont_array[i, 2] <- end_pos
                      next
                    }
                    if ((hist_cont[i] - hist_cont[i - 1]) == 
                      1) {
                      this_start <- i
                      remain_len <- cont_len - this_start
                      for (j in 1:remain_len) {
                        if ((i + j) > cont_len) {
                          break
                        }
                        if (hist_cont[i + j] - hist_cont[i + 
                          (j - 1)] == 1) {
                          end_pos <- i + j
                          cont_array[this_start:i, 2] <- end_pos
                          next
                        }
                        if ((hist_cont[i + j] - hist_cont[i + 
                          j - 1]) > 1) {
                          end_pos <- i + j
                          cont_array[this_start:i, 2] <- end_pos
                          break
                        }
                      }
                      i <- (i + j)
                    }
                    this_start <- NULL
                    start_pos <- NULL
                    end_pos <- NULL
                  }
                  cont_array[1, 2] <- cont_array[1, 1]
                  cont_array[cont_len, 2] <- cont_array[cont_len, 
                    1]
                }
                if (cont_len == 1) {
                  max_ext <- max(ext_counts)
                  max_pos <- which(ext_counts == max_ext)
                  cont_values <- ext_counts[max_pos]
                  mids_values <- mid_extr[max_pos]
                }
                if (cont_len > 1) {
                  if (hist_cont[cont_len] - hist_cont[cont_len - 
                    1] > 1) {
                    cont_array <- cont_array[-cont_len, ]
                  }
                  if (is.array(cont_array)) {
                    for (k in 1:length(cont_array[, 1])) {
                      cont_array[k, 3] <- cont_array[k, 2] - 
                        cont_array[k, 1]
                    }
                  }
                  else {
                    cont_array[3] <- cont_array[1] - cont_array[2]
                  }
                  if (is.array(cont_array)) {
                    max_value <- max(cont_array[, 3])
                    max_pos <- which(cont_array[, 3] == max_value)
                    if (length(max_pos) == 1 && cont_array[max_pos, 
                      3] > 2) {
                      max_coord <- cont_array[max_pos, 1:2]
                      cont_values <- ext_counts[max_coord[1]:max_coord[2]]
                      mids_values <- mid_extr[max_coord[1]:max_coord[2]]
                      if (length(hist_cont) > 1) {
                        if (hist_cont[cont_len] - hist_cont[cont_len - 
                          1] > 1) {
                          max_coord[2] <- (cont_array[max_pos, 
                            2]) - 1
                          cont_values <- ext_counts[max_coord[1]:max_coord[2]]
                          mids_values <- mid_extr[max_coord[1]:max_coord[2]]
                        }
                      }
                    }
                    else {
                      max_ext <- max(ext_counts)
                      max_pos <- which(ext_counts == max_ext)
                      cont_values <- ext_counts[max_pos]
                      mids_values <- mid_extr[max_pos]
                    }
                  }
                  else {
                    max_ext <- max(ext_counts)
                    max_pos <- which(ext_counts == max_ext)
                    cont_values <- ext_counts[max_pos]
                    mids_values <- mid_extr[max_pos]
                  }
                }
            }
            if (length(hist_mids) == 1) {
                cont_values <- length(hist_table)
                mids_values <- mean(hist_table)
            }
            hist_array <- array(NA, c(length(cont_values), 2))
            hist_array[, 1] <- cont_values
            hist_array[, 2] <- mids_values
            colnames(hist_array) <- c("bin_counts", "bin_midpoints")
            hist_return <- list(array = hist_array, vector = hist_table)
            return(hist_return)
        }
        PDF_fit <- function(x, normal = TRUE) {
            h <- hist(x)
            xhist <- c(min(h$breaks), h$breaks)
            yhist <- c(0, h$density, 0)
            xfit <- seq(min(x), max(x), length = 40)
            if (normal == TRUE) {
                yfit <- dnorm(xfit, mean = mean(x), sd = sd(x))
            }
            if (normal == FALSE) {
                yfit <- dlogis(xfit, location = mean(x), scale = sd(x))
            }
            plot(xhist, yhist, type = "s", ylim = c(0, max(yhist, 
                yfit)), main = "Normal pdf and histogram")
            lines(xfit, yfit, col = "red2")
        }
        TSS_clustered <- var_calc_i(TSS_data_i, gene)
        PDF_object <- cluster_PDF(TSS_clustered, binwidth = binWidth)
        hist_peak <- hist_peak_find(PDF_object, count)
        segment_vector2 <- hist_peak$vector
        vec_range <- range(hist_peak$array[, 2])
        if (length(vec_range > 1)) {
            vec_min <- vec_range[1] - (binWidth/2)
            vec_max <- vec_range[2] + (binWidth/2)
        }
        else {
            vec_min <- hist_peak$array[1, 2] - (binWidth/2)
            vec_min <- hist_peak$array[1, 2] + (binWidth/2)
        }
        segment_index <- which(segment_vector2 >= vec_min & segment_vector2 <= 
            vec_max)
        segment_vector3 <- segment_vector2[segment_index]
        TSR_shape <- kurtosis(segment_vector3)
        if (is.na(TSR_shape)) {
            TSR_shape <- c("Completely_Peaked")
        }
        TSR_breadth <- max(segment_vector3) - min(segment_vector3) + 
            1
        TSR_mid <- mean(segment_vector3)
        TSR_mid <- round(TSR_mid, 0)
        TSR_range <- range(segment_vector3)
        nTags <- length(segment_vector3)
        total_tags <- length(TSS_clustered$tss)
        peak_ratio <- nTags/total_tags
        TSR_output <- list(shape = TSR_shape, breadth = TSR_breadth, 
            range = TSR_range, mid = TSR_mid, tss_vector = segment_vector3, 
            tag_count = nTags, total_count = total_tags, Ratio = peak_ratio)
        if (draw == TRUE) {
            PDF_fit(segment_vector2, normal = TRUE)
        }
        return(TSR_output)
    }
    genes_vector <- TSS_data_i[, 2]
    genesUnique <- unique(genes_vector)
    ngenes <- length(genesUnique)
    TSR_array <- array(NA, c(ngenes, 8))
    colnames(TSR_array) <- c("Shape", "TSR_Breadth", "Peak_Count", 
        "Total_Count", "Tag_Ratio", "Midpoint", "Start", "End")
    rownames(TSR_array) <- genesUnique
    for (t in 1:ngenes) {
        this_Gene <- genesUnique[t]
        print(this_Gene)
        TSR_out <- TSR_packagev2(TSS_data, this_Gene, cluster_Num, 
            binWidth, count, Draw)
        TSR_array[t, 1] <- TSR_out$shape
        TSR_array[t, 2] <- TSR_out$breadth
        TSR_array[t, 3] <- TSR_out$tag_count
        TSR_array[t, 4] <- TSR_out$total_count
        TSR_array[t, 5] <- TSR_out$Ratio
        TSR_array[t, 6] <- TSR_out$mid
        TSR_array[t, 7:8] <- TSR_out$range
    }
    if (print) {
        write.table(TSR_array, file = "TSR_output.txt", sep = "\t", 
            row.names = TRUE, col.names = TRUE)
    }
    return(TSR_array)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
