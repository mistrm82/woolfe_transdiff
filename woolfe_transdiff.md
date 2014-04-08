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
library(mogene10sttranscriptcluster.db)
library(oligo)
library(limma)
library(plyr)
library(ggplot2)
library(grid)
library(gridExtra)
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

# Convert to factors
pData(data)$celltype <- factor(pData(data)$celltype)
pData(data)$Batch <- factor(pData(data)$Batch)

# restore previous working directory
setwd(prevWorkDir)
```



## Quality control on raw data
ArrayQualityMetrics QC report for the data can be found [here](./results/report_raw/index.html)




A heatmap cluster shows that samples cluster by cell type, which is good. Batch 2 samples tend to cluster together, but they also cluster amongst induced nociceptor samples from Batch 1.

<img src="figure/image1_.png" title="plot of chunk image1 " alt="plot of chunk image1 " width="800px" style="display: block; margin: auto;" />


Boxplots, representing summaries of the signal intensity distributions of the arrays, are mostly similar in position and widths. Arrays that look most different are in the 'induced nociceptor group' and is possibly due to the two batches. 

<img src="figure/image2_.png" title="plot of chunk image2 " alt="plot of chunk image2 " width="800px" style="display: block; margin: auto;" />


## Quality check after normalization 
The data was pre-processed and normalized using RMA() (Robust Multiarray Analysis) algorithm. RMA is the standard method for performing background subtraction, normalization and averaging. With the processed data we find that signal intensity is more similarly distributed across the arrays. Although, based on the QC report Array #8 appears to be an outlier and will be removed for downstream analysis. Density estimates show no indication of major outliers as the smoothed histograms for each array have similar shapes and ranges. 


```r
# apply rma
data.rma <- rma(data)
```

```
## Background correcting
## Normalizing
## Calculating Expression
```





<img src="figure/image3.png" title="plot of chunk image3" alt="plot of chunk image3" width="800px" />


<img src="figure/image4_.png" title="plot of chunk image4 " alt="plot of chunk image4 " width="800px" />

### PCA 
Plot all pairwise combinations of the first 5 principal components for the samples. The more similar the samples are, the closer they will cluster in these plots. In the first one below cell type is coded using the same color key as above.

<img src="figure/calculate_PCA.png" title="plot of chunk calculate_PCA" alt="plot of chunk calculate_PCA" width="800px" style="display: block; margin: auto;" />


In the second plot the color coding represents the two batches (red = batch1; cyan = batch 2). A second run of induced nociceptor samples were run in the second batch. Here, we see that while the data clusters well into three cell populations in the first and second principal components - the batches are evident in PC3. 

<img src="figure/calculate_PCA2.png" title="plot of chunk calculate_PCA2" alt="plot of chunk calculate_PCA2" width="800px" style="display: block; margin: auto;" />


## Differential expression analysis
We applied [limma](http://www.bioconductor.org/packages/release/bioc/html/limma.html) to the data and used a linear modeling approach to identify gene expression changes between the different cell populations. Both cell type and Batch were included in the model and we looked at contrasts between the different cell types. 


```r

# Load the control annotations
data(mogene10stCONTROL, package = "ArrayTools")

# Remove control probes from data
gene <- data.rma[!(row.names(exprs(data.rma)) %in% mogene10stCONTROL[, 1]), 
    ]
gene <- gene[, which(colnames(gene) %in% names(bo.out) == F)]

# Model matrix
pd <- pData(gene)
mod <- model.matrix(~celltype + Batch - 1, data = pd)

# Fit a linear model
fit <- lmFit(gene, mod)

# Define a comtrast model matrix
cont.matrix <- makeContrasts(MEFvsDRG = "celltypeMEF-celltypeDRG", MEFvsiNoc = "celltypeMEF-celltypeiNoc", 
    iNocvsDRG = "celltypeiNoc-celltypeDRG", levels = mod)

# Extract the linear model fit for the contrasts
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

# Set threshold
p.cutoff <- 0.001
logfc.cutoff <- 2

# Get gene tables
MEFvsDRG <- topTable(fit2, coef = 1, number = nrow(exprs(gene)))
MEFvsiNoc <- topTable(fit2, coef = 2, number = nrow(exprs(gene)))
iNocvsDRG <- topTable(fit2, coef = 3, number = nrow(exprs(gene)))
```


### Volcano plots
Gene expression changes were assessed by three different comparisons:

1. mouse embryonic fibrolasts vs dorsal root ganglion
2. mouse embryonic fibroblasts vs induced nociceptor 
3. induce nociceptor vs dorsal root ganglion

The third comparison is where would anticipate to see the least amount of change, under the assumption that the transdifferentiation process was successful. Along the same vein, we would probabaly expect to see a good amount of overlap between the genes identified in the first two comparisons. The volcano plots below highlight the findings from each comparison with significant probes (FDR < 0.001; logFC > 2) represented in turquoise.

<img src="figure/volcanoplot.png" title="plot of chunk volcanoplot" alt="plot of chunk volcanoplot" width="800px" style="display: block; margin: auto;" />


### Table of significant genes
<table>
 <thead>
  <tr>
   <th>  </th>
   <th> SigGenes </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td> MEFvsDRG </td>
   <td> 2885 </td>
  </tr>
  <tr>
   <td> MEFvsiNoc </td>
   <td> 1209 </td>
  </tr>
  <tr>
   <td> iNocvsDRG </td>
   <td> 1648 </td>
  </tr>
</tbody>
</table>



### Overlapping significant genes
The number of overlapping genes between the three different comparisons:
<img src="figure/image5_.png" title="plot of chunk image5 " alt="plot of chunk image5 " width="800px" />









