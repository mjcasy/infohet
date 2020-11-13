context("Het")
library(infohet)


test_that("getHet works at extremes", {
  Counts <- Matrix::sparseMatrix(i = c(1,2,2),
                                 j = c(1,1,2),
                                 x = c(10,5,5))
  expect_equal(getHet(Counts), c(1,0))
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

  Counts <- Matrix::sparseMatrix(i = c(1,1,1,1,2,2,2),
                                 j = c(1,2,3,4,1,2,3),
                                 x = c(2,2,2,2,1,1,1))
  expect_equal(getHetMicro(Counts, factor(c("1", "1", "1", "2"))), c(0,0))

})

test_that("there is one label per cell", {
  Counts <- Matrix::sparseMatrix(i = c(1,1,1,1,2,2),
                                 j = c(1,2,3,4,1,2),
                                 x = c(2,2,2,2,1,1))
  expect_error(getHetMicro(Counts, factor(c("1", "1", "2", "2", "2"))))
  expect_error(getHetMacro(Counts, factor(c("1", "1", "2", "2", "2"))))
  expect_error(groupCounts(Counts, factor(c("1", "1", "2", "2", "2"))))
})

test_that("netInformation works at extremes", {
  Counts <- Matrix::sparseMatrix(i = c(1,1,1,1,2,2),
                                 j = c(1,2,3,4,1,2),
                                 x = c(2,2,2,2,4,4))
  expect_equal(netInformation(Counts, factor(c("1", "2", "3", "4"))), 0)
  expect_equal(netInformation(Counts, factor(c("1", "1", "1", "1"))), 0)

})
