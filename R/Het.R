get_Het <- function(CountsMatrix) {
  require(Matrix)

  Total <- Matrix::rowSums(CountsMatrix)
  N <- ncol(CountsMatrix)

  CountsMatrix <- t(CountsMatrix)

  Indices <- length(CountsMatrix@p)-1
  Het <- vector("numeric", length(Indices))

  for (ind in 1:Indices) {
    count <- CountsMatrix@x[(CountsMatrix@p[ind]+1) : CountsMatrix@p[ind+1]]
    freq <- count / Total[ind]
    Het[ind] <- t(freq) %*% log(N*freq)
  }

  Het[is.infinite(Het)] <- NA
  Het[is.na(Het)] <- 0

  names(Het) <- colnames(CountsMatrix)

  Het
}

adjust_Het <- function(Het, CountsMatrix) {
  require(Matrix)

  Total <- Matrix::rowSums(CountsMatrix)
  N <- ncol(CountsMatrix)

  M <- Total
  M[M > N] <- N

  HetAdj <- Het - log(N / M)

  HetAdj
}

get_HetMacro <- function(CountsMatrix, groups, GroupedCounts) {
  require(Matrix)

  Total <- Matrix::rowSums(CountsMatrix)
  N <- ncol(CountsMatrix)

  Ng <- as.vector(table(groups))

  GroupedCounts <- t(GroupedCounts)

  Indices <- length(GroupedCounts@p)-1
  Het <- vector("numeric", length(Indices))

  for (ind in 1:Indices) {
    count <- GroupedCounts@x[(GroupedCounts@p[ind]+1) : GroupedCounts@p[ind+1]]
    freq <- count / Total[ind]

    NonZero <- 1 + GroupedCounts@i[(GroupedCounts@p[ind]+1) : GroupedCounts@p[ind+1]]

    Het[ind] <- t(freq) %*% log(N*freq /  Ng[NonZero])
  }

  Het[is.infinite(Het)] <- 0

  names(Het) <- rownames(CountsMatrix)

  Het
}

get_HetMicro <- function(CountsMatrix, groups, GroupedCounts) {
  require(Matrix)

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

  HetMicro <- cbind(HetMicro, Average)

  rownames(HetMicro) <- rownames(CountsMatrix)

  HetMicro
}
