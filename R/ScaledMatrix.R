#' The ScaledMatrix class
#'
#' Defines the ScaledMatrixSeed and ScaledMatrix classes and their associated methods.
#' These classes support delayed centering and scaling of the columns in the same manner as \code{\link{scale}},
#' but preserving the original data structure for more efficient operations like matrix multiplication.
#' 
#' @param x A matrix or any matrix-like object (e.g., from the \pkg{Matrix} package).
#' 
#' This can alternatively be a ScaledMatrixSeed, in which case any values of \code{center} and \code{scale} are ignored.
#' @param center A numeric vector of length equal to \code{ncol(x)}, where each element is to be subtracted from the corresponding column of \code{x}.
#' A \code{NULL} value indicates that no subtraction is to be performed.
#' Alternatively \code{TRUE}, in which case it is set to the column means of \code{x}.
#' @param scale A numeric vector of length equal to \code{ncol(x)}, where each element is to divided from the corresponding column of \code{x} (after subtraction).
#' A \code{NULL} value indicates that no division is to be performed.
#' Alternatively \code{TRUE}, in which case it is set to the column-wise root-mean-squared differences from \code{center} 
#' (interpretable as standard deviations if \code{center} is set to the column means, see \code{\link{scale}} for commentary).
#' 
#' @return 
#' The \code{ScaledMatrixSeed} constructor will return a ScaledMatrixSeed object.
#' 
#' The \code{ScaledMatrix} constructor will return a ScaledMatrix object equivalent to \code{t((t(x) - center)/scale)}.
#' 
#' @section Methods for ScaledMatrixSeed objects:
#' ScaledMatrixSeed objects are implemented as \linkS4class{DelayedMatrix} backends.
#' They support standard operations like \code{dim}, \code{dimnames} and \code{extract_array}.
#' 
#' Passing a ScaledMatrixSeed object to the \code{\link{DelayedArray}} constructor will create a ScaledMatrix object.
#' 
#' It is possible for \code{x} to contain a ScaledMatrix, thus nesting one ScaledMatrix inside another.
#' This can occasionally be useful in combination with transposition to achieve centering/scaling in both dimensions.
#' 
#' @section Methods for ScaledMatrix objects:
#' ScaledMatrix objects are derived from \linkS4class{DelayedMatrix} objects and support all of valid operations on the latter.
#' Several functions are specialized for greater efficiency when operating on ScaledMatrix instances, including:
#' \itemize{
#'     \item Subsetting, transposition and replacement of row/column names.
#'         These will return a new ScaledMatrix rather than a DelayedMatrix.
#'     \item Matrix multiplication via \code{\%*\%}, \code{crossprod} and \code{tcrossprod}.
#'         These functions will return a DelayedMatrix.
#'     \item Calculation of row and column sums and means by \code{colSums}, \code{rowSums}, etc. 
#' }
#' 
#' All other operations applied to a ScaledMatrix will use the underlying \pkg{DelayedArray} machinery.
#' Unary or binary operations will generally create a new DelayedMatrix instance containing a ScaledMatrixSeed.
#' 
#' Tranposition can effectively be used to allow centering/scaling on the rows if the input \code{x} is transposed.
#'
#' @section Efficiency vs precision:
#' The raison d'etre of the ScaledMatrix is that it can offer faster matrix multiplication by avoiding the \pkg{DelayedArray} block processing.
#' This is done by refactoring the scaling/centering operations to use the (hopefully more efficient) multiplication operator of the original matrix \code{x}.
#' Unfortunately, the speed-up comes at the cost of increasing the risk of catastrophic cancellation.
#' The procedure requires subtraction of one large intermediate number from another to obtain the values of the final matrix product.
#' This could result in a loss of numerical precision that compromises the accuracy of downstream algorithms. 
#' In practice, this does not seem to be a major concern though one should be careful if the input \code{x} contains very large positive/negative values.
#' 
#' @author
#' Aaron Lun
#' 
#' @examples
#' library(Matrix)
#' y <- ScaledMatrix(rsparsematrix(10, 20, 0.1), 
#'     center=rnorm(20), scale=1+runif(20))
#' y
#' 
#' crossprod(y)
#' tcrossprod(y)
#' y %*% rnorm(20)
#' 
#' @aliases
#' ScaledMatrixSeed
#' ScaledMatrixSeed-class
#' 
#' dim,ScaledMatrixSeed-method
#' dimnames,ScaledMatrixSeed-method
#' extract_array,ScaledMatrixSeed-method
#' DelayedArray,ScaledMatrixSeed-method
#' show,ScaledMatrixSeed-method
#' 
#' ScaledMatrix
#' ScaledMatrix-class
#' 
#' dimnames<-,ScaledMatrix,ANY-method
#' t,ScaledMatrix-method
#' [,ScaledMatrix,ANY,ANY,ANY-method
#' 
#' colSums,ScaledMatrix-method
#' rowSums,ScaledMatrix-method
#' colMeans,ScaledMatrix-method
#' rowMeans,ScaledMatrix-method
#' 
#' %*%,ANY,ScaledMatrix-method
#' %*%,ScaledMatrix,ANY-method
#' %*%,ScaledMatrix,ScaledMatrix-method
#' 
#' crossprod,ScaledMatrix,missing-method
#' crossprod,ScaledMatrix,ANY-method
#' crossprod,ANY,ScaledMatrix-method
#' crossprod,ScaledMatrix,ScaledMatrix-method
#' 
#' tcrossprod,ScaledMatrix,missing-method
#' tcrossprod,ScaledMatrix,ANY-method
#' tcrossprod,ANY,ScaledMatrix-method
#' tcrossprod,ScaledMatrix,ScaledMatrix-method
#' 
#' @docType class
#' @name ScaledMatrix
NULL

#' @export
#' @rdname ScaledMatrix
#' @importFrom DelayedArray DelayedArray
#' @importFrom Matrix colMeans rowSums t
ScaledMatrix <- function(x, center=NULL, scale=NULL) {
    if (isTRUE(center)) {
        center <- colMeans(x)
    }
    if (isTRUE(scale)) {
        tx <- t(DelayedArray(x))
        if (!is.null(center)) {
            tx <- tx - center
        }
        ss <- rowSums(tx^2)
        scale <- sqrt(ss / (nrow(x) - 1))
    }
    DelayedArray(ScaledMatrixSeed(x, center=center, scale=scale))
}

#' @export
#' @importFrom DelayedArray DelayedArray new_DelayedArray
setMethod("DelayedArray", "ScaledMatrixSeed",
    function(seed) new_DelayedArray(seed, Class="ScaledMatrix")
)

###################################
# Overridden utilities from DelayedArray, for efficiency.

#' @export
#' @importFrom DelayedArray DelayedArray seed
setReplaceMethod("dimnames", "ScaledMatrix", function(x, value) {
    DelayedArray(rename_ScaledMatrixSeed(seed(x), value))
})

#' @export
#' @importFrom DelayedArray DelayedArray seed
setMethod("t", "ScaledMatrix", function(x) {
    DelayedArray(transpose_ScaledMatrixSeed(seed(x)))
})

#' @export
#' @importFrom DelayedArray DelayedArray seed
setMethod("[", "ScaledMatrix", function(x, i, j, ..., drop=TRUE) {
    if (missing(i)) i <- NULL
    if (missing(j)) j <- NULL
    out <- DelayedArray(subset_ScaledMatrixSeed(seed(x), i=i, j=j))

    if (drop && any(dim(out)==1L)) {
        return(drop(out))
    }
    out
})

###################################
# Basic matrix stats.

#' @export
#' @importFrom Matrix colSums rowSums drop
setMethod("colSums", "ScaledMatrix", function(x, na.rm = FALSE, dims = 1L) {
    if (is_transposed(seed(x))) {
        return(rowSums(t(x)))
    }

    out <- rep(1, nrow(x)) %*% x
    out <- drop(out)
    names(out) <- colnames(x)
    out
})

#' @export
#' @importFrom Matrix colSums rowSums drop
setMethod("rowSums", "ScaledMatrix", function(x, na.rm = FALSE, dims = 1L) {
    if (is_transposed(seed(x))) {
        return(colSums(t(x)))
    }

    out <- x %*% rep(1, ncol(x))
    out <- drop(out)
    names(out) <- rownames(x)
    out
})

#' @export
#' @importFrom Matrix colMeans colSums
setMethod("colMeans", "ScaledMatrix", function(x, na.rm = FALSE, dims = 1L) colSums(x)/nrow(x))

#' @export
#' @importFrom Matrix rowMeans rowSums
setMethod("rowMeans", "ScaledMatrix", function(x, na.rm = FALSE, dims = 1L) rowSums(x)/ncol(x))
