context("Het")
library(infohet)



test_that("groupCounts works", {
  Counts <- Matrix::sparseMatrix(i = c(1,1,1,1,2,2),
                                 j = c(1,2,3,4,1,2),
                                 x = c(2,2,2,2,4,4))
  GroupedCounts <- Matrix::sparseMatrix(i = c(1,1,2),
                                        j = c(1,2,1),
                                        x = c(4,4,8))
  colnames(GroupedCounts) <- c("1", "2")
  expect_equal(groupCounts(Counts, factor(c("1", "1", "2", "2"))), GroupedCounts)

})


