# Tests the ScaledMatrix implementation.
# library(testthat); library(ScaledMatrix); source("setup.R"); source("test-mult.R")

##########################
# Defining a class that can't do anything but get multiplied.
# This checks that there isn't any hidden DelayedArray realization 
# happening, which would give the same results but slower.

setClass("CrippledMatrix", slots=c(x="matrix"))

setMethod("dim", c("CrippledMatrix"), function(x) dim(x@x))

setMethod("colSums", c("CrippledMatrix"), function(x) colSums(x@x))

setMethod("rowSums", c("CrippledMatrix"), function(x) rowSums(x@x))

setMethod("sweep", c("CrippledMatrix"), function (x, MARGIN, STATS, FUN = "-", check.margin = TRUE, ...) {
    sweep(x@x, MARGIN, STATS, FUN, check.margin, ...)
})

setMethod("%*%", c("CrippledMatrix", "ANY"), function(x, y) x@x %*% y)

setMethod("%*%", c("ANY", "CrippledMatrix"), function(x, y) x %*% y@x)

setMethod("crossprod", c("CrippledMatrix", "missing"), function(x, y) crossprod(x@x))

setMethod("crossprod", c("CrippledMatrix", "ANY"), function(x, y) crossprod(x@x, y))

setMethod("crossprod", c("ANY", "CrippledMatrix"), function(x, y) crossprod(x, y@x))

setMethod("tcrossprod", c("CrippledMatrix", "missing"), function(x, y) tcrossprod(x@x))

setMethod("tcrossprod", c("CrippledMatrix", "ANY"), function(x, y) tcrossprod(x@x, y))

setMethod("tcrossprod", c("ANY", "CrippledMatrix"), function(x, y) tcrossprod(x, y@x))

spawn_extra_scenarios <- function(NR=50, NC=20) {
    c(
        spawn_scenarios(NR, NC),
        spawn_scenarios_basic(NR, NC, 
            CREATOR=function(r, c) { 
                new("CrippledMatrix", x=matrix(runif(NR*NC), ncol=NC))
            }, 
            REALIZER=function(x) x@x
        )
    )
}

##########################

test_that("ScaledMatrix right multiplication works as expected", {
    possibles <- spawn_extra_scenarios(100, 50)
    for (test in possibles) {
        ref.y <- test$ref
        bs.y <- test$def

        # Multiply by a vector.
        z <- rnorm(ncol(ref.y))
        expect_equal_product(bs.y %*% z, ref.y %*% z)

        # Multiply by a matrix.
        z <- matrix(rnorm(ncol(ref.y)*10), ncol=10)
        expect_equal_product(bs.y %*% z, ref.y %*% z)

        # Multiply by an empty matrix.
        z <- matrix(0, ncol=0, nrow=ncol(ref.y))
        expect_equal_product(bs.y %*% z, ref.y %*% z)
    }
})

test_that("ScaledMatrix left multiplication works as expected", {
    possibles <- spawn_extra_scenarios(50, 80)
    for (test in possibles) {
        ref.y <- test$ref
        bs.y <- test$def

        # Multiply by a vector.
        z <- rnorm(nrow(ref.y))
        expect_equal_product(z %*% bs.y, z %*% ref.y)

        # Multiply by a matrix.
        z <- matrix(rnorm(nrow(ref.y)*10), nrow=10)
        expect_equal_product(z %*% bs.y, z %*% ref.y)

        # Multiply by an empty matrix.
        z <- matrix(0, nrow=0, ncol=nrow(ref.y))
        expect_equal_product(z %*% bs.y, z %*% ref.y)
    }
})

test_that("ScaledMatrix dual multiplication works as expected", {
    # Not using the CrippledMatrix here; some scaling of the inner matrix is unavoidable
    # when the inner matrix is _not_ a ScaledMatrix but is being multiplied by one.
    possibles1 <- spawn_scenarios(10, 20)
    for (test1 in possibles1) {
        possibles2 <- spawn_scenarios(20, 15)
        for (test2 in possibles2) {

            expect_equal_product(test1$def %*% test2$def, test1$ref %*% test2$ref)

            # Checking that zero-dimension behaviour is as expected.
            expect_equal_product(test1$def[0,] %*% test2$def, test1$ref[0,] %*% test2$ref)
            expect_equal_product(test1$def %*% test2$def[,0], test1$ref %*% test2$ref[,0])
            expect_equal_product(test1$def[,0] %*% test2$def[0,], test1$ref[,0] %*% test2$ref[0,])
            expect_equal_product(test1$def[0,] %*% test2$def[,0], test1$ref[0,] %*% test2$ref[,0])
        }
    }
})

##########################

test_that("ScaledMatrix lonely crossproduct works as expected", {
    possibles <- spawn_extra_scenarios(90, 30)
    for (test in possibles) {
        ref.y <- test$ref
        bs.y <- test$def
        expect_equal_product(crossprod(bs.y), crossprod(ref.y))
    }
})

test_that("ScaledMatrix crossproduct from right works as expected", {
    possibles <- spawn_extra_scenarios(60, 50)
    for (test in possibles) {
        ref.y <- test$ref
        bs.y <- test$def

        # Multiply by a vector.
        z <- rnorm(nrow(ref.y))
        expect_equal_product(crossprod(bs.y, z), crossprod(ref.y, z))

        # Multiply by a matrix.
        z <- matrix(rnorm(nrow(ref.y)*10), ncol=10)
        expect_equal_product(crossprod(bs.y, z), crossprod(ref.y, z))

        # Multiply by an empty matrix.
        z <- matrix(0, ncol=0, nrow=nrow(ref.y))
        expect_equal_product(crossprod(bs.y, z), crossprod(ref.y, z))
    }
})

test_that("ScaledMatrix crossproduct from left works as expected", {
    possibles <- spawn_extra_scenarios(40, 100)
    for (test in possibles) {
        ref.y <- test$ref
        bs.y <- test$def

        # Multiply by a vector.
        z <- rnorm(nrow(ref.y))
        expect_equal_product(crossprod(z, bs.y), crossprod(z, ref.y))

        # Multiply by a matrix.
        z <- matrix(rnorm(nrow(ref.y)*10), ncol=10)
        expect_equal_product(crossprod(z, bs.y), crossprod(z, ref.y))

        # Multiply by an empty matrix.
        z <- matrix(0, ncol=0, nrow=nrow(ref.y))
        expect_equal_product(crossprod(z, bs.y), crossprod(z, ref.y))
    }
})

test_that("ScaledMatrix dual crossprod works as expected", {
    possibles1 <- spawn_scenarios(20, 50)
    for (test1 in possibles1) {
        possibles2 <- spawn_scenarios(20, 15)
        for (test2 in possibles2) {

            expect_equal_product(crossprod(test1$def, test2$def), crossprod(test1$ref, test2$ref))

            # Checking that zero-dimension behaviour is as expected.
            expect_equal_product(crossprod(test1$def[,0], test2$def), crossprod(test1$ref[,0], test2$ref))
            expect_equal_product(crossprod(test1$def, test2$def[,0]), crossprod(test1$ref, test2$ref[,0]))
            expect_equal_product(crossprod(test1$def[0,], test2$def[0,]), crossprod(test1$ref[0,], test2$ref[0,]))
        }
    }
})

##########################

test_that("ScaledMatrix lonely tcrossproduct works as expected", {
    possibles <- spawn_extra_scenarios(50, 80)
    for (test in possibles) {
        ref.y <- test$ref
        bs.y <- test$def
        expect_equal_product(tcrossprod(bs.y), tcrossprod(ref.y))
    }
})

test_that("ScaledMatrix tcrossproduct from right works as expected", {
    possibles <- spawn_extra_scenarios(60, 70)
    for (test in possibles) {
        ref.y <- test$ref
        bs.y <- test$def

        # Multiply by a vector (this doesn't work).
        z <- rnorm(ncol(ref.y))
        expect_error(tcrossprod(bs.y, z), "non-conformable")
        expect_error(tcrossprod(ref.y, z), "non-conformable")

        # Multiply by a matrix.
        z <- matrix(rnorm(ncol(ref.y)*10), nrow=10)
        expect_equal_product(tcrossprod(bs.y, z), tcrossprod(ref.y, z))

        # Multiply by an empty matrix.
        z <- matrix(0, nrow=0, ncol=ncol(ref.y))
        expect_equal_product(tcrossprod(bs.y, z), tcrossprod(ref.y, z))
    }
})

test_that("ScaledMatrix tcrossproduct from left works as expected", {
    possibles <- spawn_extra_scenarios(80, 50)
    for (test in possibles) {
        ref.y <- test$ref
        bs.y <- test$def

        # Multiply by a vector.
        z <- rnorm(ncol(ref.y))
        expect_equal_product(tcrossprod(z, bs.y), tcrossprod(z, ref.y))

        # Multiply by a matrix.
        z <- matrix(rnorm(ncol(ref.y)*10), nrow=10)
        expect_equal_product(tcrossprod(z, bs.y), tcrossprod(z, ref.y))

        # Multiply by an empty matrix.
        z <- matrix(0, nrow=0, ncol=ncol(ref.y))
        expect_equal_product(tcrossprod(z, bs.y), tcrossprod(z, ref.y))
    }
})

test_that("ScaledMatrix dual tcrossprod works as expected", {
    possibles1 <- spawn_scenarios(20, 50)
    for (test1 in possibles1) {
        possibles2 <- spawn_scenarios(25, 50)
        for (test2 in possibles2) {

            expect_equal_product(tcrossprod(test1$def, test2$def), tcrossprod(test1$ref, test2$ref))

            # Checking that zero-dimension behaviour is as expected.
            expect_equal_product(tcrossprod(test1$def[0,], test2$def), tcrossprod(test1$ref[0,], test2$ref))
            expect_equal_product(tcrossprod(test1$def, test2$def[0,]), tcrossprod(test1$ref, test2$ref[0,]))
            expect_equal_product(tcrossprod(test1$def[,0], test2$def[,0]), tcrossprod(test1$ref[,0], test2$ref[,0]))
        }
    }
})

##########################

wrap_in_ScMat <- function(input, reference) 
# Wrapping an input matrix in a ScaledMatrix.
{
    output <- vector("list", 8)
    counter <- 1L

    for (trans in c(FALSE, TRUE)) {
        for (it in 1:4) {
            if (trans) { 
                y <- t(input)
                ref <- t(reference)
            } else {
                ref <- reference
                y <- input
            }

            adjusted <- scale_and_center(y, ref, it)
            if (trans) {
                adjusted$def <- t(adjusted$def)
                adjusted$ref <- t(adjusted$ref)
            }

            output[[counter]] <- adjusted
            counter <- counter+1L
        }
    }
    output
}

test_that("nested ScaledMatrix works as expected", {
    basic <- matrix(rnorm(400), ncol=20)

    available <- list(list(def=basic, ref=basic))
    for (nesting in 1:2) {
        # Creating nested ScMats with and without scaling/centering/transposition.
        next_available <- vector("list", length(available))
        for (i in seq_along(available)) {
            current <- available[[i]]
            next_available[[i]] <- wrap_in_ScMat(current$def, current$ref)
        }

        # Testing each one of the newly created ScMats.
        available <- unlist(next_available, recursive=FALSE)
        for (i in seq_along(available)) {
            test <- available[[i]]

            # Coercion works.
            expect_equal(as.matrix(test$def), test$ref)

            # Basic stats work.
            expect_equal(rowSums(test$ref), rowSums(test$def))
            expect_equal(colSums(test$ref), colSums(test$def))

            # Multiplication works.
            y <- matrix(rnorm(20*2), ncol=2)
            expect_equal_product(test$def %*% y, test$ref %*% y)
            expect_equal_product(t(y) %*% test$def, t(y) %*% test$ref)

            # Cross product.
            y <- matrix(rnorm(20*2), ncol=2)
            expect_equal_product(crossprod(test$def), crossprod(test$ref))
            expect_equal_product(crossprod(test$def, y), crossprod(test$ref, y))
            expect_equal_product(crossprod(y, test$def), crossprod(y, test$ref))

            # Transposed cross product.
            y <- matrix(rnorm(20*2), nrow=2) 
            expect_equal_product(tcrossprod(test$def), tcrossprod(test$ref))
            expect_equal_product(tcrossprod(test$def, y), tcrossprod(test$ref, y))
            expect_equal_product(tcrossprod(y, test$def), tcrossprod(y, test$ref))
        }
    }
})

set.seed(1200001)
test_that("deep testing of tcrossproduct internals: special mult", {
    NR <- 20
    NC <- 10
    basic <- matrix(rnorm(NC*NR), ncol=NC)
    c <- runif(NC)
    s <- runif(NR)

    ref <- t(matrix(c, NR, NC, byrow=TRUE)) %*% (basic/s^2)
    out <- BiocSingular:::.internal_mult_special(c, s, basic)
    expect_equal(ref, out)

    available <- list(list(def=basic, ref=basic))
    for (nesting in 1:2) {
        # Creating nested ScMats with and without scaling/centering/transposition.
        next_available <- vector("list", length(available))
        for (i in seq_along(available)) {
            current <- available[[i]]
            next_available[[i]] <- wrap_in_ScMat(current$def, current$ref)
        }

        # Testing each one of the newly created nested ScMats.
        available <- unlist(next_available, recursive=FALSE)
        for (i in seq_along(available)) {
            test <- available[[i]]
            ref <- t(matrix(c, NR, NC, byrow=TRUE)) %*% (test$ref/s^2)
            out <- BiocSingular:::.internal_mult_special(c, s, test$def)
            expect_equal(ref, out)
        }
    }
})

set.seed(1200002)
test_that("deep testing of tcrossproduct internals: scaled tcrossprod", {
    NC <- 30
    NR <- 15
    s <- runif(NC) 
    basic <- matrix(rnorm(NC*NR), ncol=NC)

    ref <- crossprod(t(basic)/s)
    out <- BiocSingular:::.internal_tcrossprod(basic, s)
    expect_equal(ref, out)

    available <- list(list(def=basic, ref=basic))
    for (nesting in 1:2) {
        # Creating nested ScMats with and without scaling/centering/transposition.
        next_available <- vector("list", length(available))
        for (i in seq_along(available)) {
            current <- available[[i]]
            next_available[[i]] <- wrap_in_ScMat(current$def, current$ref)
        }

        # Testing each one of the newly created nested ScMats.
        available <- unlist(next_available, recursive=FALSE)
        for (i in seq_along(available)) {
            test <- available[[i]]
            ref <- crossprod(t(test$ref)/s)
            out <- BiocSingular:::.internal_tcrossprod(test$def, s)
            expect_equal(ref, out)
        }
    }
})
