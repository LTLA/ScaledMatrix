# This tests the auto-scaling in the constructor.
# library(testthat); library(ScaledMatrix); source("test-scale.R")

stripscale <- function(mat, ...) {
    out <- scale(mat, ...)
    attr(out, "scaled:center") <- NULL
    attr(out, "scaled:scale") <- NULL
    out
}

test_that("ScaledMatrix mimics scale()", {
    mat <- matrix(rnorm(10000), ncol=10)

    out <- ScaledMatrix(mat, center=TRUE)
    expect_identical(as.matrix(out), stripscale(mat, scale=FALSE))

    out <- ScaledMatrix(mat, scale=TRUE)
    expect_identical(as.matrix(out), stripscale(mat, center=FALSE))

    out <- ScaledMatrix(mat, center=TRUE, scale=TRUE)
    expect_identical(as.matrix(out), stripscale(mat))

    out <- ScaledMatrix(mat)
    expect_identical(as.matrix(out), mat)
})
