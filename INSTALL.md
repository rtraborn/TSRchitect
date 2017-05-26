# __TSRchitect__ Installation and Setup

## Preliminary Steps
__TSRchitect__ is package for the [R software environment for statistical
computing](https://www.r-project.org/).
We assume that you have a current version of [R](https://www.r-project.org/)s
 installed and are familiar with its use.
Several other packages are prerequisite for __TSRchitect__, and these packages
are available through [_Bioconductor_](http://bioconductor.org/).
We assume that you have a current version of [R](https://www.r-project.org/)
installed and are familiar with its use.
Typical preliminary steps to install or updated these packages are as follows
(note: for system-wide installation, you would need to invoke
[R](https://www.r-project.org/) as root):

```{r eval=FALSE}
#installing CRAN packages
install.packages(c("gtools","knitr"))
```

```{r eval=FALSE}
source("http://bioconductor.org/biocLite.R")
biocLite(c("AnnotationHub", "BiocGenerics", "BiocParallel", "ENCODExplorer", "GenomicAlignments", "GenomeInfoDb", "GenomicRanges", "IRanges", "Rsamtools", "rtracklayer", "S4Vectors", "SummarizedExperiment"))
```

## Obtaining TSRchitect
As of April 24th, 2017, __TSRchitect__ is available as a
[_Bioconductor_](http://bioconductor.org/) package as part of the Bioconductor 3.5 release.
Package information for __TSRchitect__ on the Bioconductor site is found [here](https://www.bioconductor.org/packages/release/bioc/html/TSRchitect.html).

__TSRchitect__ can be installed directly from Bioconductor as follows:
```{r eval=TRUE}
source("http://bioconductor.org/biocLite.R")
biocLite("TSRchitect")
```

Alternatively, for the latest development version, you can install
__TSRchitect__ as follows from our github repository (again, use root privileges on the install command for
system-wide installation):

This repository contains two branches. The branch 'master' is our stable, submitted Bioconductor version, whereas our 'devel' branch contains recent development and additional files that are not found in the branch master.

If you wish to clone and install the current stable version of the software, please type the following command:
```{bash eval=FALSE}
git clone https://github.com/BrendelGroup/TSRchitect
R CMD INSTALL TSRchitect
```

To clone and install just the 'devel' branch, please do the following:
```{bash eval=FALSE}
git clone https://github.com/BrendelGroup/TSRchitect --branch devel --single-branch
R CMD INSTALL TSRchitect
```

To check on successful installation, load the package in your
R console:

```{r eval=FALSE}
library(TSRchitect)
```
Either fix problems according to displayed error messages, or go on to the next
step.

## Finally

Proceed to the vignette (inst/doc/TSRchitect.md), and we'll tell you how to execute
sample workflows (or, equally easy, your very own data analyses).