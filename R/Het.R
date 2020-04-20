#' Information in Heterogeneity
#'
#' Calculates Het, the information encoded in heterogeneity, in a gene-wise manner
#'
#' @param CountsMatrix Feature x cell sparse counts matrix of class dgCMatrix
#'
#' @return
#'
#' @export
#'
#' @details
#' Output has units of bits. Theoretical maximum of calculation is log number of cells.
#'
#' @examples
#' Counts <- Matrix::sparseMatrix(i = c(2,3,3),
#'                                j = c(1,1,2),
#'                                x = c(1,1,1))
#' Het <-  get_Het(CountsMatrix = Counts)
#'
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
#' @param Het Numeric vector of Het values of length equal to number of features
#' @param CountsMatrix Feature x cell sparse counts matrix of class dgCMatrix
#'
#' @return
#' Numeric vector of length equal to the number of features (rows in CountsMatrix)
#'
#' @export
#'
#' @details
#' Sparisty is defined here as a feature having some number of counts M, less than
#' N, the number of cells.
#' The minimum information of a sparse feature is log N/M bits, which is
#' subtracted from Het.
#'
#' @examples
#' Counts <- Matrix::sparseMatrix(i = c(1,2,2,2),
#'                                j = c(1,1,2,2),
#'                                x = rep(1, 4))
#' Het <- get_Het(Counts)
#' subtract_HetSparse(Het, Counts)
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
#' Calculates HetMacro, the heterogeneity explained by a given model of population structure,
#' in a gene-wise manner
#'
#' @param CountsMatrix Feature x cell sparse counts matrix of class dgCMatrix
#' @param groups Factor of cell identities
#' @param GroupedCounts Feature x grouped cells sparse counts matrix of class dgCMatrix
#'
#' @return
#' @export
#'
#' @examples
#' Counts <- Matrix::sparseMatrix(i = c(1,1,1,1,2,2,2),
#'                                j = c(1,2,3,4,1,2,3),
#'                                x = c(2,2,2,2,3,3,2))
#' Ident <- factor(c("1", "1", "2", "2"))
#' GroupCounts <- group_Counts(Counts, Ident)
#' get_HetMacro(Counts, Ident, GroupCounts)
get_HetMacro <- function(CountsMatrix, groups, GroupedCounts) {

  if(length(groups) != ncol(CountsMatrix)){
    warning("Inconsistent number of cells between objects:\n\tlength(groups) != ncol(CountsMatrix)")
  }

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
#' Calculates HetMicro, the heterogeneity within each group of cells, or its weigthed
#' average (default)
#'
#' @param CountsMatrix Feature x cell sparse counts matrix of class dgCMatrix
#' @param groups Factor of cell identities
#' @param GroupedCounts Feature x grouped cells sparse counts matrix of class dgCMatrix
#' @param full Logical flag for whether to return just HetMicro or full Het by group
#'
#' @return
#' @export
#'
#' @examples
#' Counts <- Matrix::sparseMatrix(i = c(1,1,1,1,2,2,2),
#'                                j = c(1,2,3,4,1,2,3),
#'                                x = c(2,2,2,2,3,3,2))
#' Ident <- factor(c("1", "1", "2", "2"))
#' GroupCounts <- group_Counts(Counts, Ident)
#' get_HetMicro(Counts, Ident, GroupCounts)
get_HetMicro <- function(CountsMatrix, groups, GroupedCounts, full = F) {

  if(length(groups) != ncol(CountsMatrix)){
    warning("Inconsistent number of cells between objects:\n\tlength(groups) != ncol(CountsMatrix)")
  }

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

  if(full == F){
    HetMicro <- Average
    names(HetMicro) <- rownames(CountsMatrix)
  } else if(full == T){
    HetMicro <- cbind(HetMicro, Average)
    rownames(HetMicro) <- rownames(CountsMatrix)
  }

  HetMicro
}
