
<!-- README.md is generated from README.Rmd. Please edit that file -->

# infohet

<!-- badges: start -->

<!-- badges: end -->

Package for the quantification of the information content of single-cell RNA-sequencing data-sets. Based on the quantification of information in heterogeneity (infohet).

Below are example applications to feature selection, cluster assessment and identification of differentially expressed genes.

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
library(ggplot2)

# Sparse Counts Matrix
load("../Data/10x/CountsMatrix")

# Factor of cell identities (i.e. cluster labels)
load("../Data/10x/Identity")

minTotal <- 100
infoThreshold <- 0.5

Total <- Matrix::rowSums(CountsMatrix)
if(any(Total < minTotal)){
   CountsMatrix <- CountsMatrix[-which(Total < minTotal),]
   Total <- Total[-which(Total < minTotal)]
}
```

Feature Selection. Select genes with the most information (Het) unexplained by technical effects (Hom) or sparsity (HetSparse). Technical model samples from categorical distribution in proportion to cell library depths. 

``` r
Het <- getHet(CountsMatrix)
HetAdj <- subtractHetSparse(Het, CountsMatrix)

nullHet <- simulateHom(CountsMatrix)
nullHet <- subtractHetSparse(nullHet, CountsMatrix)

Threshold <- nullHet+infoThreshold

N <- CountsMatrix@Dim[2]
Mean_nUMI <- Total / N

HetDataFrame <- data.frame(log10(Mean_nUMI), HetAdj, nullHet, Threshold, HetAdj > Threshold)
colnames(HetDataFrame) <- c("log10_Mean_nUMI", "Het", "Null_Model", "Threshold", "Selected")

ggplot(HetDataFrame, aes(x = log10_Mean_nUMI, y = Het, colour = Selected)) + geom_point() +
  geom_line(aes(y = Null_Model), colour = "black") + 
  ylim(0, log2(N))
#> Warning: Removed 14 rows containing missing values (geom_point).
#> Warning: Removed 484 row(s) containing missing values (geom_path).
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

Cluster Assessment. A good clustering of cells should minimise unexplained information (HetMicro).

``` r
GroupedCounts <- groupCounts(CountsMatrix, Identity)

HetMicro <- getHetMicro(CountsMatrix, Identity, GroupedCounts)
HetMicroAdj <- subtractHetSparse(HetMicro, CountsMatrix)

HetDataFrame <- cbind(HetDataFrame, HetMicroAdj)

ggplot(HetDataFrame, aes(x = log10_Mean_nUMI, y = HetMicroAdj, colour = Selected)) + geom_point() +
  geom_line(aes(y = Null_Model), colour = "black") + 
  ylim(-log2(minTotal), log2(N))
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

Differential Expression. Identification of genes with substantial information explainable by one cluster vs remaining (HetMacro). Most genes will have a non-zero amount of HetMacro due to technical effects. Either an arbitrary threshold (e.g. 0.05 bits) should be used or the labelling can be permuted and the mean HetMacro from a set of permutations used.

``` r
DEGroup <- factor(ifelse(Identity == 17, "1", "2"))

GroupedCounts <- groupCounts(CountsMatrix, DEGroup)

HetMacro <- getHetMacro(CountsMatrix, DEGroup, GroupedCounts)

HetDataFrame <- cbind(HetDataFrame, HetMacro)

ggplot(HetDataFrame, aes(x = log10_Mean_nUMI, y = HetMacro, colour = HetMacro > 0.05)) + geom_point() +
  ylim(0, log2(N))
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" />
