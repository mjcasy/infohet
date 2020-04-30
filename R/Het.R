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
#' Het <-  getHet(CountsMatrix = Counts)
#'
getHet <- function(CountsMatrix) {

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
#' Sparsity is defined here as a feature having some number of counts M, less than
#' N, the number of cells.
#' The minimum information of a sparse feature is log N/M bits, which is
#' subtracted from Het.
#'
#' @examples
#' Counts <- Matrix::sparseMatrix(i = c(1,2,2,2),
#'                                j = c(1,1,2,2),
#'                                x = rep(1, 4))
#' Het <- getHet(Counts)
#' subtractHetSparse(Het, Counts)
subtractHetSparse <- function(Het, CountsMatrix) {

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
#' @param Groups Factor of cell identities
#' @param groupedCounts Feature x grouped cells sparse counts matrix of class dgCMatrix
#'
#' @return
#' @export
#'
#' @examples
#' Counts <- Matrix::sparseMatrix(i = c(1,1,1,1,2,2,2),
#'                                j = c(1,2,3,4,1,2,3),
#'                                x = c(2,2,2,2,3,3,2))
#' Ident <- factor(c("1", "1", "2", "2"))
#' GroupCounts <- groupCounts(Counts, Ident)
#' getHetMacro(Counts, Ident, GroupCounts)
getHetMacro <- function(CountsMatrix, Groups, groupedCounts) {

  if(length(Groups) != ncol(CountsMatrix)){
    warning("Inconsistent number of cells between objects:\n\tlength(Groups) != ncol(CountsMatrix)")
  }

  Total <- Matrix::rowSums(CountsMatrix)
  N <- ncol(CountsMatrix)

  Ng <- as.vector(table(Groups))

  groupedCounts <- Matrix::t(groupedCounts)

  Indices <- length(groupedCounts@p)-1
  Het <- vector("numeric", length(Indices))

  for (ind in 1:Indices) {
    count <- groupedCounts@x[(groupedCounts@p[ind]+1) : groupedCounts@p[ind+1]]
    freq <- count / Total[ind]

    NonZero <- 1 + groupedCounts@i[(groupedCounts@p[ind]+1) : groupedCounts@p[ind+1]]

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
#' @param Groups Factor of cell identities
#' @param groupedCounts Feature x grouped cells sparse counts matrix of class dgCMatrix
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
#' GroupCounts <- groupCounts(Counts, Ident)
#' getHetMicro(Counts, Ident, GroupCounts)
getHetMicro <- function(CountsMatrix, Groups, groupedCounts, full = F) {

  if(length(Groups) != ncol(CountsMatrix)){
    warning("Inconsistent number of cells between objects:\n\tlength(Groups) != ncol(CountsMatrix)")
  }

  Total <- Matrix::rowSums(CountsMatrix)
  N <- ncol(CountsMatrix)

  types <- levels(Groups)
  I <- length(types)

  CountsList <- list()

  for (i in 1:I) {
    CountsList[[i]] <- CountsMatrix[,Groups == types[i], drop = FALSE]
  }

  HetList <- lapply(CountsList, getHet)

  HetMicro <- matrix(unlist(HetList), nrow = nrow(CountsMatrix), ncol = I)
  HetMicro[is.infinite(HetMicro)] <- 0
  colnames(HetMicro) <- types
  rownames(HetMicro) <- rownames(CountsMatrix)

  FreqMatrix <- groupedCounts * Total^-1

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
