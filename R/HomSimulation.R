
#' Simulate Homogenous gene expression
#'
#' Simulates null model of gene expression for given count matrix, producing
#' values for Het expected for each gene under the null.
#'
#' @param CountsMatrix Feature x cell sparse counts matrix of class dgCMatrix
#' @param NumTrials Number of trials of simulation
#' @param SimStep Step size, in log10 counts, between number of counts simulated
#' @param Spline Logical flag for whether to simulate each gene directly or infer Het from splining
#' @param DepthAdjusted Logical flag for whether simulated cells should be
#'                       assigned counts in ratio of corresponding library sizes
#' @param subtractSparsity Subtract information due to count sparsity
#'
#' @return Numeric vector of gene-wise simulated Het
#' @export
#'
#' @details
#' Null model is simulated over set of cells of size N, where N is the number of cells in CountsMatrix.
#' Some number of counts is randomly assigned to cells either with equal probability or in proportion to observed count depths (total counts per cell)
#' Het is calculated for each trial. The mean Het for each number of counts is found.
#' The null Het for each gene is then assigned by interpolation.
#'
#'
#' @examples
#' Counts <- Matrix::sparseMatrix(i = c(1,1,2,2),
#'                                j = c(1,2,1,2),
#'                                x = c(1000,9000,5000,5000))
#' simulateHom(Counts)
simulateHom <- function(CountsMatrix, NumTrials = 50, SimStep = 0.2, Spline = T, DepthAdjusted = T, subtractSparsity = F) {

  Total <- Matrix::rowSums(CountsMatrix)
  Range <- range(log10(Total))
  if(Spline == F){
    TotalCounts <- Total
  } else if(Spline == T){
    TotalCounts <- round(10^seq(1, Range[2], SimStep))
  }

  CellTotal <- Matrix::colSums(CountsMatrix)
  Probs <- CellTotal / sum(CellTotal)

  TrialDepths <- rep(TotalCounts, each = NumTrials)
  NumCells <- ncol(CountsMatrix)

  if (DepthAdjusted == T) {
    SimHet <- sapply(TrialDepths, function(NumTranscripts){
      Cells <- sample(NumCells, NumTranscripts, replace = T, prob = Probs)
      Ones <- rep(1, NumTranscripts)
      Mat <- Matrix::sparseMatrix(j = Cells, i = Ones, x = Ones)
      getHet(Mat)
    })
  } else {
    SimHet <- sapply(TrialDepths, function(NumTranscripts){
      Cells <- sample(NumCells, NumTranscripts, replace = T)
      Ones <- rep(1, NumTranscripts)
      Mat <- Matrix::sparseMatrix(j = Cells, i = Ones, x = Ones)
      getHet(Mat)
    })
  }

  meanHet <- tapply(SimHet, as.factor(TrialDepths), mean)

  if(Spline == F){
    NullHet <- meanHet[as.character(Total)]
  } else if(Spline == T){
    Nullfun <- stats::splinefun(TotalCounts, meanHet)
    NullHet <- Nullfun(Total)
  }

  names(NullHet) <- rownames(CountsMatrix)

  if(subtractSparsity == T){
    NullHet <- subtractHetSparse(CountsMatrix, NullHet)
  }

  NullHet
}
