# TSRchitect
## Identifying promoters from large-scale 5' end data in plants and other eukaryotes

### Installation
#### First, packages from CRAN
> install.packages(c("moments","stats"))
> library(c('moments','stats'))

#### Now packages from Bioconductor
> source("http://bioconductor.org/biocLite.R")
> biocLite(c("GenomicRanges","biomaRt","Sushi","GenomicFeatures"))
> biocLite("TSRchitect")

### Load TSRchitect
> library(TSRchitect)

### Checking out the vignette for TSRchitect
> vignette("TSRchitect")