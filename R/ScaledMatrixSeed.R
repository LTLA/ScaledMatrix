#' @export
#' @importFrom methods new is
ScaledMatrixSeed <- function(x, center=NULL, scale=NULL) {
    if (missing(x)) {
        x <- matrix(0, 0, 0)
    } else if (is(x, "ScaledMatrixSeed")) {
        return(x)
    } 

    use_center <- !is.null(center)
    use_scale <- !is.null(scale)
    new("ScaledMatrixSeed", .matrix=x, center=as.numeric(center), scale=as.numeric(scale), use_center=use_center, use_scale=use_scale, transposed=FALSE)
}

#' @importFrom S4Vectors setValidity2
setValidity2("ScaledMatrixSeed", function(object) {
    msg <- character(0)

    # Checking scalars.
    if (length(use_center(object))!=1L) {
        msg <- c(msg, "'use_center' must be a logical scalar")
    } 
    if (length(use_scale(object))!=1L) {
        msg <- c(msg, "'use_scale' must be a logical scalar")
    } 
    if (length(is_transposed(object))!=1L) {
        msg <- c(msg, "'transposed' must be a logical scalar")
    } 

    # Checking vectors.
    if (use_center(object) && length(get_center(object))!=ncol(object)) {
        msg <- c(msg, "length of 'center' must equal 'ncol(object)'")
    }
    if (use_scale(object) && length(get_scale(object))!=ncol(object)) {
        msg <- c(msg, "length of 'scale' must equal 'ncol(object)'")
    }

    if (length(msg)) {
        return(msg)
    } 
    return(TRUE)
})

#' @export
#' @importFrom methods show
setMethod("show", "ScaledMatrixSeed", function(object) {
    cat(sprintf("%i x %i ScaledMatrixSeed object", nrow(object), ncol(object)),
        sprintf("representation: %s", class(get_matrix2(object))),
        sprintf("centering: %s", if (use_center(object)) "yes" else "no"),
        sprintf("scaling: %s", if (use_scale(object)) "yes" else "no"),
    sep="\n")
})

###################################
# Internal getters. 

get_matrix2 <- function(x) x@.matrix

get_center <- function(x) x@center

get_scale <- function(x) x@scale

use_center <- function(x) x@use_center

use_scale <- function(x) x@use_scale

is_transposed <- function(x) x@transposed

###################################
# DelayedArray support utilities. 

#' @export
setMethod("dim", "ScaledMatrixSeed", function(x) {
    d <- dim(get_matrix2(x))
    if (is_transposed(x)) { d <- rev(d) }
    d
})

#' @export
setMethod("dimnames", "ScaledMatrixSeed", function(x) {
    d <- dimnames(get_matrix2(x))
    if (is_transposed(x)) { d <- rev(d) }
    d
})

#' @export
#' @importFrom DelayedArray extract_array
setMethod("extract_array", "ScaledMatrixSeed", function(x, index) {
    x2 <- subset_ScaledMatrixSeed(x, index[[1]], index[[2]])        
    realize_ScaledMatrixSeed(x2)
})

###################################
# Other utilities. 

rename_ScaledMatrixSeed <- function(x, value) {
    if (is_transposed(x)) value <- rev(value)
    dimnames(x@.matrix) <- value
    x
}

transpose_ScaledMatrixSeed <- function(x) {
    x@transposed <- !is_transposed(x)
    x
}

#' @importFrom Matrix t as.matrix
#' @importFrom methods is
realize_ScaledMatrixSeed <- function(x, ...) {
    out <- get_matrix2(x)

    if (use_scale(x) || use_center(x)) {
        if (is(out, "ScaledMatrix")) {
            # Any '-' and '/' would collapse this to a DelayedArray, 
            # which would then call extract_array, which would then 
            # call realize_ScaledMatrixSeed, forming an infinite loop.
            # So we might as well realize it now.
            out <- realize_ScaledMatrixSeed(seed(out))
        }

        out <- t(out)
        if (use_center(x)) {
            out <- out - get_center(x)
        }
        if (use_scale(x)) {
            out <- out / get_scale(x)
        }

        if (!is_transposed(x)) out <- t(out)
    } else {
        if (is_transposed(x)) out <- t(out) 
    }

    as.matrix(out)
}

subset_ScaledMatrixSeed <- function(x, i, j) {
    if (is_transposed(x)) {
        x2 <- transpose_ScaledMatrixSeed(x)
        x2 <- subset_ScaledMatrixSeed(x2, i=j, j=i)
        return(transpose_ScaledMatrixSeed(x2))
    }

    if (!is.null(i)) {
        x@.matrix <- get_matrix2(x)[i,,drop=FALSE]
    }
    
    if (!is.null(j)) {
        if (is.character(j)) {
            j <- match(j, colnames(x))
        }
        
        x@.matrix <- get_matrix2(x)[,j,drop=FALSE]
            
        if (use_scale(x)) {
            x@scale <- get_scale(x)[j]
        }
        
        if (use_center(x)) {
            x@center <- get_center(x)[j]
        }
    }

    return(x)
}
