#' Information in Heterogeneity
#'
#' Calculates Het, the information encoded in heterogeneity, in a gene-wise manner
#'
#' @param CountsMatrix Feature x cell sparse counts matrix of class dgCMatrix
#' @param base base of log
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
getHet <- function(CountsMatrix, base = 2) {

  Total <- Matrix::rowSums(CountsMatrix)
  N <- ncol(CountsMatrix)

  transposeCounts <- Matrix::t(CountsMatrix)

  Indices <- length(transposeCounts@p)-1
  Het <- vector("numeric", length(Indices))

  for (ind in 1:Indices) {
    freqshrink <- getFreqShrink(transposeCounts, ind, N, Total)
    Log2Nfreq <- log(N*freqshrink, base = base)
    Log2Nfreq[Log2Nfreq == -Inf] <- 0
    Het[ind] <- t(freqshrink) %*% Log2Nfreq
  }


  Het[is.infinite(Het)] <- 0

  names(Het) <- colnames(transposeCounts)

  Het
}

#' Calaculate Information in Model
#'
#' Calculates HetMacro, the heterogeneity explained by a given model of population structure,
#' in a gene-wise manner
#'
#' @param CountsMatrix Feature x cell sparse counts matrix of class dgCMatrix
#' @param Groups Factor of cell identities
#' @param base base of log
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
getHetMacro <- function(CountsMatrix, Groups, base = 2) {

  if(length(Groups) != ncol(CountsMatrix)){
    stop("Inconsistent number of cells between objects:\n\tlength(Groups) != ncol(CountsMatrix)")
  }

  Total <- Matrix::rowSums(CountsMatrix)
  N <- ncol(CountsMatrix)

  Ng <- as.vector(table(Groups))

  transposeCounts <- Matrix::t(CountsMatrix)

  Indices <- length(transposeCounts@p)-1
  Het <- vector("numeric", length(Indices))

  for (ind in 1:Indices) {
    freqshrink <- getFreqShrink(transposeCounts, ind, N, Total)
    groupedfreqshrink <- tapply(freqshrink, Groups, sum)

    NonZero <- which(groupedfreqshrink != 0)

    Het[ind] <- t(groupedfreqshrink[NonZero]) %*% log(N*groupedfreqshrink[NonZero] /  Ng[NonZero], base = base)
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
#' @param base base of log
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
getHetMicro <- function(CountsMatrix, Groups, base = 2) {

  if(length(Groups) != ncol(CountsMatrix)){
    stop("Inconsistent number of cells between objects:\n\tlength(Groups) != ncol(CountsMatrix)")
  }

  Total <- Matrix::rowSums(CountsMatrix)
  N <- ncol(CountsMatrix)

  types <- levels(Groups)
  I <- length(types)

  Ng <- as.vector(table(Groups))
  names(Ng) <- names(table(Groups))

  transposeCounts <- Matrix::t(CountsMatrix)

  Indices <- length(transposeCounts@p)-1
  HetMicro <- vector("numeric", length =  nrow(CountsMatrix))
  names(HetMicro) <- rownames(CountsMatrix)

  for (ind in 1:Indices) {
    freqshrink <- getFreqShrink(transposeCounts, ind, N, Total)
    groupedfreqshrink <- tapply(freqshrink, Groups, sum)
    typeHet <- vector(mode = "numeric", length = I)
    names(typeHet) <- types
    for(type in types){
      freq <- freqshrink[Groups == type] / groupedfreqshrink[type]
      typeHet[type] <- t(freq) %*% log(Ng[type]*freq, base = base)
    }
    HetMicro[ind] <- groupedfreqshrink[types] %*% typeHet
  }

  HetMicro
}

#' Calculate net total information explained
#'
#' Calculates total macro-heterogeneity regularised by the entropy of the group structure
#'
#' @param CountsMatrix Feature x cell sparse counts matrix of class dgCMatrix
#' @param Groups Factor of cell identities
#' @param Het Vector of feature Het values
#' @param lambdaMod Scale factor of regularisation; shoud be in the range 0 to 1
#' @param base base of log
#'
#' @return Net total information explained
#' @export
#'
#' @details The default regularisation assumes that most heterogeneity is going
#' to be explained by differences between clusters. If there is substantial
#' intra-cluster heterogeneity relative to inter-cluster differences, negative
#' net information may be returned even for `good' clusterings. A lower
#' lambdaMod value should be used in such cases, but not so low that the net
#' information remains high/positive for excessive numbers of clusters.
#'
#' @examples
#' Counts <- Matrix::sparseMatrix(i = c(1,1,1,1,2,2,2),
#'                                j = c(1,2,3,4,1,2,3),
#'                                x = c(2,2,2,2,3,3,2))
#' Ident <- factor(c("1", "1", "2", "2"))
#' netInformation(Counts, Ident)
netInformation <- function(CountsMatrix, Groups, Het, lambdaMod = 1, base = 2){

  if(lambdaMod < 0 | lambdaMod > 1){
    stop("lambdaMod should be between 0 and 1")
  }

  HetMacro <- getHetMacro(CountsMatrix, Groups, base = base)

  N <- ncol(CountsMatrix)
  g <- nrow(CountsMatrix)

  fS <- as.vector(table(Groups)) / N
  HS <- as.numeric(-1 * t(fS) %*% log(fS, base = base))

  if(missing(Het)){
    Het <- getHet(CountsMatrix, base = base)
  }
  stopifnot(length(Het) == length(HetMacro))
  lambda <- sum(Het) / log(N, base = base)

  nI_S <- sum(HetMacro) - lambdaMod*lambda*HS

  nI_S

}
