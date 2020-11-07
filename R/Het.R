#' Information in Heterogeneity
#'
#' Calculates Het, the information encoded in heterogeneity, in a gene-wise manner
#'
#' @param CountsMatrix Feature x cell sparse counts matrix of class dgCMatrix
#' @param shrinkage Boolean flag on whether to use James-Stein type shrinkage estimator
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
getHet <- function(CountsMatrix, shrinkage = T, subtractSparsity = F) {

  Total <- Matrix::rowSums(CountsMatrix)
  N <- ncol(CountsMatrix)

  transposeCounts <- Matrix::t(CountsMatrix)

  Indices <- length(transposeCounts@p)-1
  Het <- vector("numeric", length(Indices))

  if(shrinkage == T) {
    for (ind in 1:Indices) {
      freqshrink <- getFreqShrink(transposeCounts, ind, N, Total)
      Log2Nfreq <- log2(N*freqshrink)
      Log2Nfreq[Log2Nfreq == -Inf] <- 0
      Het[ind] <- t(freqshrink) %*% Log2Nfreq
    }
  } else if(shrinkage == F) {
    for (ind in 1:Indices) {
      count <- transposeCounts@x[(transposeCounts@p[ind]+1) : transposeCounts@p[ind+1]]
      freq <- count / Total[ind]
      Het[ind] <- t(freq) %*% log2(N*freq)
    }
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
#' @param shrinkage Boolean flag on whether to use James-Stein type shrinkage estimator
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
getHetMacro <- function(CountsMatrix, Groups, shrinkage = T) {

  if(length(Groups) != ncol(CountsMatrix)){
    stop("Inconsistent number of cells between objects:\n\tlength(Groups) != ncol(CountsMatrix)")
  }

  Total <- Matrix::rowSums(CountsMatrix)
  N <- ncol(CountsMatrix)

  Ng <- as.vector(table(Groups))

  if(shrinkage == T){
    transposeCounts <- Matrix::t(CountsMatrix)

    Indices <- length(transposeCounts@p)-1
    Het <- vector("numeric", length(Indices))

    for (ind in 1:Indices) {
      freqshrink <- getFreqShrink(transposeCounts, ind, N, Total)
      groupedfreqshrink <- tapply(freqshrink, Groups, sum)

      NonZero <- which(groupedfreqshrink != 0)

      Het[ind] <- t(groupedfreqshrink[NonZero]) %*% log2(N*groupedfreqshrink[NonZero] /  Ng[NonZero])
    }
  } else if(shrinkage == F){
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
#' @param shrinkage Boolean flag on whether to use James-Stein type shrinkage estimator
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
getHetMicro <- function(CountsMatrix, Groups, shrinkage = T, full = F, subtractSparsity = F, components = F) {

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

  if(shrinkage == T){
    Ng <- as.vector(table(Groups))
    names(Ng) <- names(table(Groups))

    transposeCounts <- Matrix::t(CountsMatrix)

    Indices <- length(transposeCounts@p)-1
    HetMicro <- matrix(nrow =  Indices, ncol = I, dimnames = list(rownames(CountsMatrix), types))
    fullMicro <- matrix(nrow =  Indices, ncol = I, dimnames = list(rownames(CountsMatrix), types))

    for (ind in 1:Indices) {
      freqshrink <- getFreqShrink(transposeCounts, ind, N, Total)
      groupedfreqshrink <- tapply(freqshrink, Groups, sum)
      typeHet <- vector(mode = "numeric", length = I)
      names(typeHet) <- types
      for(type in types){
        freq <- freqshrink[Groups == type] / groupedfreqshrink[type]
        typeHet[type] <- t(freq) %*% log2(Ng[type]*freq)
      }
      fullMicro[ind,] <- typeHet
      HetMicro[ind,] <- groupedfreqshrink[types] * typeHet
    }

    if(full == F & components == F){
      HetMicro <- rowSums(HetMicro)
    } else if(full == T){
      HetMicro <- fullMicro
    }

  } else if(shrinkage == F){
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
  }

  HetMicro
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


#' Calculate regularised total information explained
#'
#' Calculates total macro-heterogeneity regularised by the entropy of the group structure
#'
#' @param CountsMatrix Feature x cell sparse counts matrix of class dgCMatrix
#' @param Groups Factor of cell identities
#' @param shrinkage Boolean flag on whether to use James-Stein type shrinkage estimator
#' @param lambda Scale factor of regularisation; defaults to ratio of total information obtained to theoretical maximum
#'
#' @return Regularised total information explained
#' @export
#'
#' @examples
#' Counts <- Matrix::sparseMatrix(i = c(1,1,1,1,2,2,2),
#'                                j = c(1,2,3,4,1,2,3),
#'                                x = c(2,2,2,2,3,3,2))
#' Ident <- factor(c("1", "1", "2", "2"))
#' regularisedTotalHetMacro(Counts, Ident)
regularisedTotalHetMacro <- function(CountsMatrix, Groups, shrinkage = T, lambda){

  HetMacro <- getHetMacro(CountsMatrix, Groups, shrinkage = shrinkage)

  N <- ncol(CountsMatrix)
  g <- nrow(CountsMatrix)

  fS <- as.vector(table(Groups)) / N
  HS <- as.numeric(-1 * t(fS) %*% log2(fS))

  if(missing(lambda)){
    Het <- getHet(CountsMatrix, shrinkage = shrinkage)
    lambda <- sum(Het) / (g*log2(N))
  }

  rI_S <- sum(HetMacro) - lambda*g*HS

  rI_S

}
