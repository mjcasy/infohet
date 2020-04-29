
#' Simulate Homogenous gene expression
#'
#' Simulates null model of gene expression for given count matrix, producing
#' values for Het expected for each gene under the null.
#'
#' @param CountsMatrix Feature x cell sparse counts matrix of class dgCMatrix
#' @param NumTrials Number of trials of simulation
#' @param SimStep Step size, in log10 counts, between number of counts simulated
#' @param TotalCounts Optional integer vector of total counts for simualted cells.
#'                    Default set internally is 10 to max cell counts with SimStep steps
#' @param depth_adjusted Logical flag for whether simulated cells should be
#'                       assigned counts in ratio of corresponding library sizes
#'
#' @return
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
#' SimCounts <- 10000
#' simulate_Hom(Counts)
simulate_Hom <- function(CountsMatrix, NumTrials = 50, SimStep = 0.2, TotalCounts = NA, depth_adjusted = T) {

  Totals <- Matrix::rowSums(CountsMatrix)
  Range <- range(log10(Totals))
  TotalCounts <- round(10^seq(1, Range[2], SimStep))

  CellTotal <- Matrix::colSums(CountsMatrix)
  Probs <- CellTotal / sum(CellTotal)

  TrialDepths <- rep(TotalCounts, each = NumTrials)
  NumCells <- ncol(CountsMatrix)

  if (depth_adjusted == T) {
    SimHet <- sapply(TrialDepths, function(NumTranscripts){
      Cells <- sample(NumCells, NumTranscripts, replace = T, prob = Probs)
      Ones <- rep(1, NumTranscripts)
      Mat <- Matrix::sparseMatrix(j = Cells, i = Ones, x = Ones)
      get_Het(Mat)
    })
  } else {
    SimHet <- sapply(TrialDepths, function(NumTranscripts){
      Cells <- sample(NumCells, NumTranscripts, replace = T)
      Ones <- rep(1, NumTranscripts)
      Mat <- Matrix::sparseMatrix(j = Cells, i = Ones, x = Ones)
      get_Het(Mat)
    })
  }

  meanHet <- tapply(SimHet, as.factor(TrialDepths), mean)
  Nullfun <- stats::splinefun(TotalCounts, meanHet)

  NullHet <- Nullfun(Totals)
  names(NullHet) <- rownames(CountsMatrix)

  NullHet
}
