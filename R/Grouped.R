
#' Sum counts by group
#'
#' For a given grouping of cells, produces a sparse matrix where each column is the total counts in each group.
#'
#' @param CountsMatrix Feature x cell sparse counts matrix of class dgCMatrix
#' @param groups Factor of cell identities
#'
#' @return
#' @export
#'
#' @examples
#' Counts <- Matrix::sparseMatrix(i = c(1,1,1,1,2,2,2),
#'                                j = c(1,2,3,4,1,2,3),
#'                                x = c(2,2,2,2,3,3,2))
#' Ident <- factor(c("1", "1", "2", "2"))
#' group_Counts(Counts, Ident)
group_Counts <- function(CountsMatrix, groups) {

  if(length(groups) != ncol(CountsMatrix)){
    warning("Inconsistent number of cells between objects:\n\tlength(groups) != ncol(CountsMatrix)")
  }

  types <- levels(groups)
  I <- length(types)

  CountsList <- list()

  for (i in 1:I) {
    CountsList[[i]] <- CountsMatrix[,groups == types[i], drop = FALSE]
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
