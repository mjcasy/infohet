
#' Sum counts by group
#'
#' For a given grouping of cells, produces a sparse matrix where each column is the total counts for each group.
#'
#' @param CountsMatrix Feature x cell sparse counts matrix of class dgCMatrix
#' @param Groups Factor of cell identities
#'
#' @return Feature x group matrix of class dgCMatrix
#' @export
#'
#' @examples
#' Counts <- Matrix::sparseMatrix(i = c(1,1,1,1,2,2,2),
#'                                j = c(1,2,3,4,1,2,3),
#'                                x = c(2,2,2,2,3,3,2))
#' Ident <- factor(c("1", "1", "2", "2"))
#' groupCounts(Counts, Ident)
groupCounts <- function(CountsMatrix, Groups) {

  if(length(Groups) != ncol(CountsMatrix)){
    stop("Inconsistent number of cells between objects:\n\tlength(Groups) != ncol(CountsMatrix)")
  }

  types <- levels(Groups)
  I <- length(types)

  CountsList <- list()

  for (i in 1:I) {
    CountsList[[i]] <- CountsMatrix[,Groups == types[i], drop = FALSE]
  }

  GroupedCounts <- lapply(CountsList, Matrix::rowSums)

  GroupedCounts <- Matrix::Matrix(
    data = unlist(GroupedCounts),
    nrow = nrow(CountsMatrix),
    ncol = I,
    dimnames = list(dimnames(CountsMatrix)[[1]], types),
    sparse = TRUE
  )

  GroupedCounts
}

#' Calculate cell frequencies for shrinakge estimator
#'
#' @param transposeCounts Transposed sparse count matrix
#' @param ind Integer indicating chosen gene (row number in count matrix)
#' @param N Number of cells
#' @param Total Integer of total counts per cell
#'
#' @return Numeric vector of shrinkage cell frequencies
#' @export
#'
#' @examples
#' Counts <- Matrix::sparseMatrix(i = c(2,3,3),
#'                                j = c(1,1,2),
#'                                x = c(1,1,1))
#' N <- ncol(Counts)
#' Total <- Matrix::rowSums(Counts)
#' getFreqShrink(Matrix::t(Counts), 1, N, Total)
#'
getFreqShrink <- function(transposeCounts, ind, N, Total){

  tk <- rep(1/N, N)
  tkadj <- tk

  count <- transposeCounts@x[(transposeCounts@p[ind]+1) : transposeCounts@p[ind+1]]
  elements <- transposeCounts@i[(transposeCounts@p[ind]+1) : transposeCounts@p[ind+1]]+1

  freq <- count / Total[ind]

  num <- 1 - sum(freq^2)
  tkadj[elements] <- tkadj[elements] - freq
  den <- (sum(count) - 1)*sum(tkadj^2)
  lambda <- num/den
  lambda[lambda > 1] <- 1

  freqshrink <- lambda*tk
  freqshrink[elements]  <- freqshrink[elements] + (1 - lambda)*freq

  freqshrink
}
