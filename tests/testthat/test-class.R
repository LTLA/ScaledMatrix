# Tests the ScaledMatrix implementation.
# library(testthat); library(ScaledMatrix); source("setup.R"); source("test-class.R")

set.seed(100001)
test_that("ScaledMatrix utility functions work as expected", {
    possibles <- spawn_scenarios()
    for (test in possibles) {
        expect_s4_class(test$def, "ScaledMatrix")
        expect_identical(test$def, ScaledMatrix(DelayedArray::seed(test$def)))

        expect_identical(dim(test$def), dim(test$ref))
        expect_identical(extract_array(test$def, list(1:10, 1:10)), test$ref[1:10, 1:10])
        expect_identical(extract_array(test$def, list(1:10, NULL)), test$ref[1:10,])
        expect_identical(extract_array(test$def, list(NULL, 1:10)), test$ref[,1:10])
        expect_identical(as.matrix(test$def), test$ref)

        expect_equal(rowSums(test$def), rowSums(test$ref))
        expect_equal(colSums(test$def), colSums(test$ref))
        expect_equal(rowMeans(test$def), rowMeans(test$ref))
        expect_equal(colMeans(test$def), colMeans(test$ref))

        tdef <- t(test$def)
        expect_s4_class(tdef, "ScaledMatrix") # still a DefMat!
        expect_identical(t(tdef), test$def)
        expect_identical(as.matrix(tdef), t(test$ref))

        # Checking column names getting and setting.
        spawn_names <- sprintf("THING_%i", seq_len(ncol(test$def)))
        colnames(test$def) <- spawn_names
        expect_identical(spawn_names, colnames(test$def))
        expect_s4_class(test$def, "ScaledMatrix") # still a DefMat!
    }
})

set.seed(10000101)
test_that("ScaledMatrix silly inputs work as expected", {
    default <- ScaledMatrix()
    expect_identical(dim(default), c(0L, 0L))
    val <- as.matrix(default)
    dimnames(val) <- NULL
    expect_identical(val, matrix(0,0,0))

    # Checking erronious inputs.
    y <- matrix(rnorm(400), ncol=20)
    expect_error(ScaledMatrix(y, center=1), "length of 'center' must equal")
    expect_error(ScaledMatrix(y, scale=1), "length of 'scale' must equal")
})

set.seed(1000011)
test_that("ScaledMatrix subsetting works as expected", {
    expect_identical_and_defmat <- function(x, y) {
        expect_s4_class(x, "ScaledMatrix") # class is correctly preserved by direct seed modification.
        expect_identical(as.matrix(x), y)
    }

    possibles <- spawn_scenarios()
    for (test in possibles) {
        i <- sample(nrow(test$def))
        j <- sample(ncol(test$def))
        expect_identical_and_defmat(test$def[i,], test$ref[i,])
        expect_identical_and_defmat(test$def[,j], test$ref[,j])
        expect_identical_and_defmat(test$def[i,j], test$ref[i,j])
        
        # Works with zero dimensions.
        expect_identical_and_defmat(test$def[0,], test$ref[0,])
        expect_identical_and_defmat(test$def[,0], test$ref[,0])
        expect_identical_and_defmat(test$def[0,0], test$ref[0,0])
        
        # Dimension dropping works as expected.
        expect_identical(test$def[i[1],], test$ref[i[1],])
        expect_identical(test$def[,j[1]], test$ref[,j[1]])
        expect_identical_and_defmat(test$def[i[1],drop=FALSE], test$ref[i[1],,drop=FALSE])
        expect_identical_and_defmat(test$def[,j[1],drop=FALSE], test$ref[,j[1],drop=FALSE])

        # Transposition with subsetting works as expected.
        alt <- t(test$def)
        expect_identical(t(alt[,i]), test$def[i,])
        expect_identical(t(alt[j,]), test$def[,j])

        # Subsetting behaves with column names.
        spawn_names <- sprintf("THING_%i", seq_len(ncol(test$def)))
        colnames(test$def) <- spawn_names
        colnames(test$ref) <- spawn_names
        ch <- sample(spawn_names)
        expect_identical_and_defmat(test$def[,ch], test$ref[,ch])
    }
})

test_that("DelayedMatrix wrapping works", {
    possibles <- spawn_scenarios(80, 50)
    for (test in possibles) {
        expect_equal_product(test$def+1, test$ref+1)

        v <- rnorm(nrow(test$def))
        expect_equal_product(test$def+v, test$ref+v)
        expect_equal_product(test$def*v, test$ref*v)

        w <- rnorm(ncol(test$def))
        expect_equal_product(sweep(test$def, 2, w, "*"), sweep(test$ref, 2, w, "*"))
    }
})
