Woolfe Lab: Transdifferentiation of neurons
========================================================





Array analysis for Liz Buttermore (elizabeth.buttermore@childrens.harvard.edu) at Woolfe Lab. Contact Meeta Mistry (mmistry@hsph.harvard.edu) for additional details. Request from client was:
  
> Basically we are comparing expression profiles of three populations of cells and there is a good amount of variability between replicate samples for one of the populations because it is a fibroblast-derived neuronal cell
line. We would like to go a little further with the data than we have so far

## Bioconductor and R libraries used

```r

library(arrayQualityMetrics)
library(RColorBrewer)
library(pd.mogene.1.0.st.v1)
library(oligo)
library(limma)
library(plyr)
library(ggplot2)
library(grid)
library(xtable)

source("~/R/scripts/useful functions/roc.R")
```


### Get variables
- get base directory for analyses
- specify data and results directories
- specify column headers used in metadata file



```r
# Setup directory variables
baseDir <- "."
dataDir <- file.path(baseDir, "data")
metaDir <- file.path(dataDir, "meta")
resultsDir <- file.path(baseDir, "results")
```


### Load the data and metadata

```r
prevWorkDir = getwd()

setwd("data/")
pd <- read.AnnotatedDataFrame("covars.desc")
data <- read.celfiles(rownames(pData(pd)))
```

```
## Reading in : DRG 1__28MoGene-1_0-st-v1_29.CEL
## Reading in : DRG 2 _28MoGene-1_0-st-v1_29.CEL
## Reading in : DRG 3 _28MoGene-1_0-st-v1_29.CEL
## Reading in : MEF 1 _28MoGene-1_0-st-v1_29.CEL
## Reading in : MEF 2 _28MoGene-1_0-st-v1_29.CEL
## Reading in : MEF 3 _28MoGene-1_0-st-v1_29.CEL
## Reading in : Mouse iNoc 2_2B _28MoGene-1_0-st-v1_29.CEL
## Reading in : Mouse iNoc 3_2B _28MoGene-1_0-st-v1_29.CEL
## Reading in : Mouse iNoc 4_2B _28MoGene-1_0-st-v1_29.CEL
## Reading in : Mouse iNoc 1_2B_batch2 _8MoGene-1_0-st-v1_29.CEL
## Reading in : Mouse iNoc 2_2B_batch2 _28MoGene-1_0-st-v1_29.CEL
## Reading in : Mouse iNoc 3_2B_batch2 _28MoGene-1_0-st-v1_29.CEL
```

```r
phenoData(data) <- pd

# restore previous working directory
setwd(prevWorkDir)
```



## Quality control on raw data
ArrayQualityMetrics QC report for the data can be found [here](./results/report_raw/index.html)




from the QC report, the data look ok. A heatmap cluster shows that samples cluster by cell type, whch is good. Batch 2 samples are somewhat close to one another, but its possible this will be less evident after normalization.

<img src="figure/image1_.png" title="plot of chunk image1 " alt="plot of chunk image1 " width="800px" style="display: block; margin: auto;" />


Boxplots representing summaries of the signal intensity distributions of the arrays, are also mostly similar in position in widths. Arrays that look most different are in the 'induced nociceptor group' and is possibly due to the two batches.

<img src="figure/image2_.png" title="plot of chunk image2 " alt="plot of chunk image2 " width="800px" style="display: block; margin: auto;" />


## Quality check after normalization 
The data was pre-processed and normalized using RMA (Robust Multiarray Analysis) algorithm. RMA is the standard method for performing background subtraction, normalization and averaging. With the processed data we find that signal intensity is more simialrly distributed across the arrays. Density estimates also show no indication of lareg outliers as the smoothed histograms for each array have similar shapes and ranges. 



```r

# apply rma
data.rma <- rma(data)
```

```
## Background correcting
## Normalizing
## Calculating Expression
```

```r

# plot some quick quality figures
boxplot(data.rma, main = "", xaxt = "n")
```

<img src="figure/qcnorm1.png" title="plot of chunk qcnorm" alt="plot of chunk qcnorm" width="800px" />

```r
hist(data.rma, main = "")
```

<img src="figure/qcnorm2.png" title="plot of chunk qcnorm" alt="plot of chunk qcnorm" width="800px" />






