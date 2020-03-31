
#' Sum counts by group
#'
#' For a given grouping of cells, produces a sparse matrix where each column is the total counts in each group.
#'
#' @param CountsMatrix dgCMatrix
#' @param groups factor
#'
#' @return dgCMatrix
#' @export
#'
#' @examples
#'
group_Counts <- function(CountsMatrix, groups) {

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
