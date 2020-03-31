
#' Simualtion of Homogenous gene expression
#'
#' @param CountsMatrix dgCMatrix
#' @param TotalCounts vector
#' @param NumTrials integer
#' @param depth_adjusted logical
#'
#' @return
#' @export
#'
#' @examples
simulate_Hom <- function(CountsMatrix, TotalCounts, NumTrials, depth_adjusted = T) {

  Totals <- Matrix::rowSums(CountsMatrix)

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
