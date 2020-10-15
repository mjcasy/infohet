#' Information in Heterogeneity
#'
#' Calculates Het, the information encoded in heterogeneity, in a gene-wise manner
#'
#' @param CountsMatrix Feature x cell sparse counts matrix of class dgCMatrix
#' @param subtractSparsity Subtract information due to count sparsity
#'
#' @return Numeric vector of gene-wise Het
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
getHet <- function(CountsMatrix, subtractSparsity = F) {

  Total <- Matrix::rowSums(CountsMatrix)
  N <- ncol(CountsMatrix)

  transposeCounts <- Matrix::t(CountsMatrix)

  Indices <- length(transposeCounts@p)-1
  Het <- vector("numeric", length(Indices))

  for (ind in 1:Indices) {
    count <- transposeCounts@x[(transposeCounts@p[ind]+1) : transposeCounts@p[ind+1]]
    freq <- count / Total[ind]
    Het[ind] <- t(freq) %*% log2(N*freq)
  }

  Het[is.infinite(Het)] <- 0

  names(Het) <- colnames(transposeCounts)

  if(subtractSparsity == T){
    Het <- subtractHetSparse(CountsMatrix, Het)
  }

  Het
}

#' Subtract information due to count sparsity
#'
#' @param CountsMatrix Feature x cell sparse counts matrix of class dgCMatrix
#' @param Het Numeric vector of Het values of length equal to number of features
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
#' subtractHetSparse(Counts, Het)
subtractHetSparse <- function(CountsMatrix, Het) {

  Total <- Matrix::rowSums(CountsMatrix)
  N <- ncol(CountsMatrix)

  M <- Total
  M[M > N] <- N

  HetAdj <- Het - log2(N / M)

  HetAdj[M == 0] <- NA

  HetAdj
}

#' Calaculate Information in Model
#'
#' Calculates HetMacro, the heterogeneity explained by a given model of population structure,
#' in a gene-wise manner
#'
#' @param CountsMatrix Feature x cell sparse counts matrix of class dgCMatrix
#' @param Groups Factor of cell identities
#'
#' @return Numeric vector of gene-wise HetMacro
#' @export
#'
#' @examples
#' Counts <- Matrix::sparseMatrix(i = c(1,1,1,1,2,2,2),
#'                                j = c(1,2,3,4,1,2,3),
#'                                x = c(2,2,2,2,3,3,2))
#' Ident <- factor(c("1", "1", "2", "2"))
#' getHetMacro(Counts, Ident)
getHetMacro <- function(CountsMatrix, Groups) {

  if(length(Groups) != ncol(CountsMatrix)){
    stop("Inconsistent number of cells between objects:\n\tlength(Groups) != ncol(CountsMatrix)")
  }

  Total <- Matrix::rowSums(CountsMatrix)
  N <- ncol(CountsMatrix)

  Ng <- as.vector(table(Groups))

  groupedCounts <- groupCounts(CountsMatrix, Groups)

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
#' @param full Logical flag for whether to return just HetMicro or full Het by group
#' @param subtractSparsity Subtract information due to count sparsity. If full also TRUE, also applies to each group.
#' @param components Logical flag for whether to return just HetMicro or additive components of HetMicro by group
#'
#' @return Numeric vector of gene-wise HetMicro
#' @export
#'
#' @examples
#' Counts <- Matrix::sparseMatrix(i = c(1,1,1,1,2,2,2),
#'                                j = c(1,2,3,4,1,2,3),
#'                                x = c(2,2,2,2,3,3,2))
#' Ident <- factor(c("1", "1", "2", "2"))
#' getHetMicro(Counts, Ident)
getHetMicro <- function(CountsMatrix, Groups, full = F, subtractSparsity = F, components = F) {

  if(length(Groups) != ncol(CountsMatrix)){
    stop("Inconsistent number of cells between objects:\n\tlength(Groups) != ncol(CountsMatrix)")
  }

  if(sum(full, components) == 2){
    stop("Both components and full should not be TRUE")
  }

  Total <- Matrix::rowSums(CountsMatrix)
  N <- ncol(CountsMatrix)

  types <- levels(Groups)
  I <- length(types)

  CountsList <- list()

  for (i in 1:I) {
    CountsList[[i]] <- CountsMatrix[,Groups == types[i], drop = FALSE]
  }

  HetList <- lapply(X = CountsList, FUN = getHet, subtractSparsity = F)

  HetMicro <- matrix(unlist(HetList), nrow = nrow(CountsMatrix), ncol = I)
  HetMicro[is.na(HetMicro)] <- 0
  HetMicro[is.infinite(HetMicro)] <- 0
  colnames(HetMicro) <- types
  rownames(HetMicro) <- rownames(CountsMatrix)

  groupedCounts <- groupCounts(CountsMatrix, Groups)

  FreqMatrix <- groupedCounts * Total^-1

  Components <- FreqMatrix * HetMicro

  Overall <- Matrix::rowSums(Components)

  if(subtractSparsity == T){
    for(i in 1:ncol(HetMicro)){
      HetMicro[,i] <- subtractHetSparse(CountsMatrix = CountsList[[i]], Het = HetMicro[,i])
    }
    Overall <- subtractHetSparse(CountsMatrix, Overall)
  }

  if(full == F & components == F){
    HetMicro <- Overall
    names(HetMicro) <- rownames(CountsMatrix)
  } else if(full == T){
    HetMicro <- cbind(HetMicro, Overall)
    rownames(HetMicro) <- rownames(CountsMatrix)
  } else if(components == T){
    HetMicro <- cbind(as.matrix(Components), Overall)
    rownames(HetMicro) <- rownames(CountsMatrix)
  }

  HetMicro
}
