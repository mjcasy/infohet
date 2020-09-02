context("Het")
library(infohet)


test_that("getHet works at extremes", {
  Counts <- Matrix::sparseMatrix(i = c(2,3,3),
                                 j = c(1,1,2),
                                 x = c(1,1,1))
  expect_equal(getHet(Counts), c(0,1,0))
  expect_equal(getHet(Counts, subtractSparsity = T), c(-Inf,0,0))
})

test_that("getHetMacro works at extremes", {
  Counts <- Matrix::sparseMatrix(i = c(1,1,1,1,2,2),
                                 j = c(1,2,3,4,1,2),
                                 x = c(2,2,2,2,4,4))
  expect_equal(getHetMacro(Counts, factor(c("1", "1", "2", "2"))), getHet(Counts))

  Counts <- Matrix::sparseMatrix(i = c(1,1,1,1,2,2,2),
                                 j = c(1,2,3,4,1,2,3),
                                 x = c(2,2,2,2,4,4,4))
  expect_equal(getHetMacro(Counts, factor(c("1", "1", "1", "2"))), getHet(Counts))

})

test_that("getHetMicro works at extremes", {
  Counts <- Matrix::sparseMatrix(i = c(1,1,1,1,2,2),
                                 j = c(1,2,3,4,1,2),
                                 x = c(2,2,2,2,1,1))
  expect_equal(getHetMicro(Counts, factor(c("1", "1", "2", "2"))), c(0,0))

  Counts <- Matrix::sparseMatrix(i = c(1,1,1,1,2,2),
                                 j = c(1,2,3,4,1,2),
                                 x = c(2,2,2,2,1,1))
  expect_equal(getHetMicro(Counts, factor(c("1", "1", "2", "2")), subtractSparsity = T), c(0,-1))

  Counts <- Matrix::sparseMatrix(i = c(1,1,1,1,2,2,2),
                                 j = c(1,2,3,4,1,2,3),
                                 x = c(2,2,2,2,1,1,1))
  expect_equal(getHetMicro(Counts, factor(c("1", "1", "1", "2"))), c(0,0))

  Counts <- Matrix::sparseMatrix(i = c(1,1,1,1,2,2,2),
                                 j = c(1,2,3,4,1,2,3),
                                 x = c(2,2,2,2,1,1,1))
  expect_equal(getHetMicro(Counts, factor(c("1", "1", "1", "2")), subtractSparsity = T), c(0,-log2(4/3)))

})

test_that("there is one label per cell", {
  Counts <- Matrix::sparseMatrix(i = c(1,1,1,1,2,2),
                                 j = c(1,2,3,4,1,2),
                                 x = c(2,2,2,2,1,1))
  expect_error(getHetMicro(Counts, factor(c("1", "1", "2", "2", "2"))))
  expect_error(getHetMacro(Counts, factor(c("1", "1", "2", "2", "2"))))
  expect_error(groupCounts(Counts, factor(c("1", "1", "2", "2", "2"))))
})

test_that("full and components flags work", {
  Counts <- Matrix::sparseMatrix(i = c(1,1,1,1,2,2),
                                 j = c(1,2,3,4,1,3),
                                 x = c(2,2,2,2,2,4))

  cTarget <- c(1/3, 2/3)
  names(cTarget) <- c("1", "2")
  fTarget <- c(1, 1)
  names(fTarget) <- c("1", "2")

  expect_equal(getHetMicro(Counts, factor(c("1", "1", "2", "2")), components = T)[2,1:2], cTarget)
  expect_equal(getHetMicro(Counts, factor(c("1", "1", "2", "2")), full = T)[2,1:2], fTarget)

  expect_error(getHetMicro(Counts, factor(c("1", "1", "2", "2")), full = T, components = T))

})
