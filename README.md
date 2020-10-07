
<!-- README.md is generated from README.Rmd. Please edit that file -->

# infohet

<!-- badges: start -->

[![Codecov test
coverage](https://codecov.io/gh/mjcasy/infohet/branch/master/graph/badge.svg)](https://codecov.io/gh/mjcasy/infohet?branch=master)
<!-- badges: end -->

Package for the quantification of the information content of single-cell
RNA-sequencing data-sets, and how much of this information has been
explained by clustering. Based on the quantification of information in
heterogeneity (infohet). Preprint -
<https://www.biorxiv.org/content/10.1101/2020.10.01.322255v1>

## Installation

``` r
install.packages("devtools")
devtools::install_github("mjcasy/infohet")
```

## Workflow

General setup. Load in Data and filter low expressing genes (less than
100 transcripts total).

``` r
library(infohet)
library(RColorBrewer)

load("../Data/Tian2018/CountsMatrix")

infoThreshold <- 0.5
minTotal <- 100

Total <- Matrix::rowSums(CountsMatrix)
if(any(Total < minTotal)){
   CountsMatrix <- CountsMatrix[-which(Total < minTotal),]
   Total <- Total[-which(Total < minTotal)]
}
```

Visualise gene-wise information - the amount of information gain from
knowing the cellular allocation of transcripts for each gene. This
corresponds to the information left unexplained by assuming each gene is
homogenous in expression.

Homogeneity, adjusted for the difference in total count depths of cells,
is simulated to provide a null baseline of information.

``` r
Het <- getHet(CountsMatrix)

nullHet <- simulateHom(CountsMatrix)

HighlyInformative <- Het > nullHet + infoThreshold

N <- CountsMatrix@Dim[2]
Mean_nUMI <- Total / N

HetDataFrame <- data.frame(log10(Mean_nUMI), Het, nullHet)
colnames(HetDataFrame) <- c("log10_Mean_nUMI", "Unexplained_Information", "Null_Model")

ColourSelected <- brewer.pal(9, "Blues")[8]
ColourNotSelected <- brewer.pal(9, "Blues")[4]
Order <- order(HetDataFrame$log10_Mean_nUMI)

with(HetDataFrame, plot(log10_Mean_nUMI, 
                        Unexplained_Information,
                        col = ifelse(HighlyInformative, ColourSelected, ColourNotSelected),
                        ylim = c(0, log2(N)),
                        pch = 20)
     )
lines(HetDataFrame$log10_Mean_nUMI[Order], HetDataFrame$Null_Model[Order], col = "red", lwd = 2)
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

The example Tian dataset has a known cluster structure. The amount of
information left unexplained by assuming each gene is homogenously
expressed within each cluster is found.

``` r
load("../Data/Tian2018/Identity")

Clustered_Unexplained_Information <- getHetMicro(CountsMatrix, Identity)
HetDataFrame <- cbind(HetDataFrame, Clustered_Unexplained_Information)

HighlyInformativeDecomp <- Clustered_Unexplained_Information > nullHet + infoThreshold

with(HetDataFrame, plot(log10_Mean_nUMI, 
                        Clustered_Unexplained_Information,
                        col = ifelse(HighlyInformativeDecomp, ColourSelected, ColourNotSelected),
                        ylim = c(0, log2(N)),
                        pch = 20)
     )
lines(HetDataFrame$log10_Mean_nUMI[Order], HetDataFrame$Null_Model[Order], col = "red", lwd = 2)
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

Information is additive between independent sources. Taking each gene as
independent, the total information left unexplained in the transcriptome
can be found.

The goal of clustering is to identify a set of homogenous
subpopulations, i.e. cell types. Therefore, the chosen clustering should
minimise the amount of information left unexplained by the assumption of
cluster homogeneity.

All else being equal, information will never decrease with increasing
cluster number. There is therefore a trade-off - the optimal clustering
is that which leaves the least information unexplained in the fewest
clusters. This can be assessed by an elbow analysis, shown below using
the Seurat (v3) clustering pipeline.

The chosen clustering can be reapplied to identify which genes are well
explained by the clustering (see above figure).

``` r
library(Seurat)

SeuObj <- Seurat::CreateSeuratObject(CountsMatrix)
SeuObj <- Seurat::SCTransform(SeuObj)

SeuObj <- RunPCA(SeuObj, verbose = FALSE)

SeuObj <- FindNeighbors(SeuObj, dims = 1:30, verbose = FALSE)

NumClusters <- c()
InformationUnexplained <- c()
Resolutions <- c(seq(0.0001, 0.001, 0.0001), seq(0.002, 0.01, 0.001), seq(0.02, 0.2, 0.01), seq(0.3, 1, 0.1))

for(i in 1:length(Resolutions)){
  Idents(SeuObj) <- NA
  SeuObj <- FindClusters(SeuObj, verbose = FALSE, resolution = Resolutions[i])
  
  Identity <- Idents(SeuObj)
  
  HetMicro <- getHetMicro(CountsMatrix, Identity)
  
  NumClusters[i] <- length(levels(Identity))
  InformationUnexplained[i] <- sum(HetMicro) 
}

Elbow <- cbind(Resolutions, NumClusters, InformationUnexplained)
```

Elbow Plot of total information explained by clustering against
resolution hyperparameter and cluster
number

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" /><img src="man/figures/README-unnamed-chunk-7-2.png" width="100%" />
