
<!-- README.md is generated from README.Rmd. Please edit that file -->

# infohet

<!-- badges: start -->

<!-- badges: end -->

There is a need for the robust quantification of cellular heterogeneity
in single-cell RNA-sequencing data. This package is based around an
information-theoretic measure of heterogeneity, quantifying
heterogeneity as the information required to produce the observed
pattern of gene expression (Het).

Cellular heterogeneity can be broadly split into that due to the
presence of multiple distinct cell types (macro-heterogeneity) and that
due to technical effects and stochastic fluctuations
(micro-heterogeneity). Information is additively decomposable so that
for a set labelling of cells, it can be split into that information
explainable by the labelling (HetMacro), and that left unexplained
(HetMicro).

\(Het(X) = HetMacro(X,l) + HetMicro(X,l)\)

Where X is the observed gene expression distribution and l is a
labelling of cells, for example by cell type.

We can apply this framework to various scRNA-seq analysis tasks,
including feature selection, cluster assessment and identification of
differentially expressed genes.

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

Feature Selection based on Het, the information content in gene
expression. The Het of each gene is found and adjusted for sparsity.
Genes with excessive Het compared to simulation of the null are
identified for selection. The null distribution is either a discrete
uniform or a multinomial with probabilities in proportion to cell count
depths (default).

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

Cluster Quality based on HetMicro. HetMicro is the gene-wise measure of
information left unexplained by some labelling of cells. This labelling
is typically the results of clustering for cell type identification.
Genes with excessive HetMicro are inadequately explained by the
clustering.

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

Differential Gene Expression based on HetMacro. HetMacro is the
information explained by a labelling of cells. When only two unique
labels are supplied, HetMacro is analogous to DGE tests. Most genes will
have a non-zero amount of HetMacro due to technical effects. Either an
arbitrary threshold (e.g. 0.05 bits) should be used or the labelling can
be permuted and the mean HetMacro from a set of permutations used.

``` r
DEGroup <- factor(ifelse(Identity == 17, "1", "2"))

GroupedCounts <- groupCounts(CountsMatrix, DEGroup)

HetMacro <- getHetMacro(CountsMatrix, DEGroup, GroupedCounts)

HetDataFrame <- cbind(HetDataFrame, HetMacro)

ggplot(HetDataFrame, aes(x = log10_Mean_nUMI, y = HetMacro, colour = HetMacro > 0.05)) + geom_point() +
  ylim(0, log2(N))
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" />
