context("HomSimulation")
library(infohet)

test_that("extremes of simulation approximately correct", {
  Counts <- Matrix::sparseMatrix(i = c(1,1,2,2),
                                 j = c(1,2,1,2),
                                 x = c(2*10^6,10,0,1))
  SetSimStep <- log10(2*10^6) - 1

  set.seed(1)
  T1 <- simulateHom(Counts, SimStep = SetSimStep)
  set.seed(1)
  T2 <- simulateHom(Counts, DepthAdjusted = F, SimStep = SetSimStep)
  expect_equal(T1[1], 1, tolerance = 1e-1)
  expect_equal(T1[2], 0, tolerance = 1e-4)
  expect_equal(T2[1], 0, tolerance = 1e-4)
  expect_true(T2[2] > 1e-4)

  expect_equal(subtractHetSparse(Counts, T1)[2], T1[2] - log2(2/1))
  expect_equal(subtractHetSparse(Counts, T2)[2], T2[2] - log2(2/1))
})

test_that("simulation sparsity correction works", {
  Counts <- Matrix::sparseMatrix(i = c(1,1,2,2),
                                 j = c(1,2,1,2),
                                 x = c(10,10,0,1))

  set.seed(1)
  T3 <- simulateHom(Counts)

  M <- Matrix::rowSums(Counts)
  M[M > 2] <- 2

  set.seed(1)
  expect_equal(simulateHom(Counts, subtractSparsity = T), T3 - log2(2/M))
})

test_that("extremes of simulation approximately correct without spline", {
  Counts <- Matrix::sparseMatrix(i = c(1,1,2,2),
                                 j = c(1,2,1,2),
                                 x = c(2*10^6,10,0,1))
  set.seed(1)
  T4 <- simulateHom(Counts, Spline = F)
  expect_equal(T4[1], 1, tolerance = 1e-1)
  expect_equal(T4[2], 0, tolerance = 1e-4)
})

