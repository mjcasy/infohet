group_Counts <- function(CountsMatrix, groups) {
  require(Matrix)

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
