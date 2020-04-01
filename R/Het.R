#' Calculate Information in Heterogeneity
#'
#' @param CountsMatrix dgCMatrix
#'
#' @return
#' @export
#'
#' @examples
get_Het <- function(CountsMatrix) {

  Total <- Matrix::rowSums(CountsMatrix)
  N <- ncol(CountsMatrix)

  CountsMatrix <- Matrix::t(CountsMatrix)

  Indices <- length(CountsMatrix@p)-1
  Het <- vector("numeric", length(Indices))

  for (ind in 1:Indices) {
    count <- CountsMatrix@x[(CountsMatrix@p[ind]+1) : CountsMatrix@p[ind+1]]
    freq <- count / Total[ind]
    Het[ind] <- t(freq) %*% log2(N*freq)
  }

  Het[is.infinite(Het)] <- NA
  Het[is.na(Het)] <- 0

  names(Het) <- colnames(CountsMatrix)

  Het
}

#' Subtract information due to count sparsity
#'
#' @param Het vector
#' @param CountsMatrix dgCMatrix
#'
#' @return
#' @export
#'
#' @examples
subtract_HetSparse <- function(Het, CountsMatrix) {

  Total <- Matrix::rowSums(CountsMatrix)
  N <- ncol(CountsMatrix)

  M <- Total
  M[M > N] <- N

  HetAdj <- Het - log2(N / M)

  HetAdj
}

#' Calaculate Information in Model
#'
#' @param CountsMatrix dgCMatrix
#' @param groups factor
#' @param GroupedCounts dgCMatrix
#'
#' @return
#' @export
#'
#' @examples
get_HetMacro <- function(CountsMatrix, groups, GroupedCounts) {

  Total <- Matrix::rowSums(CountsMatrix)
  N <- ncol(CountsMatrix)

  Ng <- as.vector(table(groups))

  GroupedCounts <- Matrix::t(GroupedCounts)

  Indices <- length(GroupedCounts@p)-1
  Het <- vector("numeric", length(Indices))

  for (ind in 1:Indices) {
    count <- GroupedCounts@x[(GroupedCounts@p[ind]+1) : GroupedCounts@p[ind+1]]
    freq <- count / Total[ind]

    NonZero <- 1 + GroupedCounts@i[(GroupedCounts@p[ind]+1) : GroupedCounts@p[ind+1]]

    Het[ind] <- t(freq) %*% log2(N*freq /  Ng[NonZero])
  }

  Het[is.infinite(Het)] <- 0

  names(Het) <- rownames(CountsMatrix)

  Het
}

#' Calaculate Information not in Model
#'
#' @param CountsMatrix dgCMatrix
#' @param groups factor
#' @param GroupedCounts dgCMatrix
#' @param reduced logical
#'
#' @return
#' @export
#'
#' @examples
get_HetMicro <- function(CountsMatrix, groups, GroupedCounts, reduced = T) {

  Total <- Matrix::rowSums(CountsMatrix)
  N <- ncol(CountsMatrix)

  types <- levels(groups)
  I <- length(types)

  CountsList <- list()

  for (i in 1:I) {
    CountsList[[i]] <- CountsMatrix[,groups == types[i], drop = FALSE]
  }

  HetList <- lapply(CountsList, get_Het)

  HetMicro <- matrix(unlist(HetList), nrow = nrow(CountsMatrix), ncol = I)
  HetMicro[is.infinite(HetMicro)] <- 0
  colnames(HetMicro) <- types
  rownames(HetMicro) <- rownames(CountsMatrix)

  FreqMatrix <- GroupedCounts * Total^-1

  Average <- Matrix::rowSums(FreqMatrix * HetMicro)

  if(reduced == T){
    HetMicro <- Average
    names(HetMicro) <- rownames(CountsMatrix)
  } else if(reduced == F){
    HetMicro <- cbind(HetMicro, Average)
    rownames(HetMicro) <- rownames(CountsMatrix)
  }

  HetMicro
}
