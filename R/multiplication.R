# We attempt to use operators defined for '.matrix' in the 'ScaledMatrixSeed'.
# This avoids expensive modifications such as loss of sparsity.
# Centering and scaling are factored out into separate operations.
#
# We assume that the non-'ScaledMatrix' argument is small and can be modified cheaply.
# We also assume that the matrix product is small and can be modified cheaply.
# This allows centering and scaling to be applied *after* multiplication.
#
# Here are some ground rules for how these functions must work:
#
#  - NO arithmetic operations shall be applied to a ScaledMatrix.
#    This includes nested ScaledMatrices that are present in '.matrix'.
#    Such operations collapses the ScaledMatrix to a DelayedMatrix, 
#    resulting in slow block processing during multiplication.
# 
#  - NO addition/subtraction operations shall be applied to '.matrix'.
#    This is necessary to avoid loss of sparsity for sparse '.matrix',
#    as well as to avoid block processing for ScaledMatrix '.matrix'.
# 
#  - NO division/multiplication operations should be applied to '.matrix'.
#    This is largely a consequence of the first point above.
#    Exceptions are only allowed when this is unavoidable, e.g., in '.internal_tcrossprod'.
#
#  - NO calling of %*% or (t)crossprod on a ScaledMatrix of the same nesting depth as an input ScaledMatrix.
#    Internal multiplication should always be applied to '.matrix', to avoid infinite S4 recursion.
#    Each method call should strip away one nesting level, i.e., operate on the seed.
#    Exceptions are allowed for dual ScaledMatrix multiplication,
#    where one argument is allowed to be of the same depth.

#' @export
#' @importFrom Matrix t
#' @importFrom DelayedArray seed DelayedArray
setMethod("%*%", c("ScaledMatrix", "ANY"), function(x, y) {
    x_seed <- seed(x)
    if (is_transposed(x_seed)) {
        out <- t(.leftmult_ScaledMatrix(t(y), x_seed))
    } else {
        out <- .rightmult_ScaledMatrix(x_seed, y)
    }
    DelayedArray(out)
})

#' @importFrom DelayedArray sweep
.rightmult_ScaledMatrix <- function(x_seed, y) {
    if (use_scale(x_seed)) {
        y <- y / get_scale(x_seed)
    }

    out <- as.matrix(get_matrix2(x_seed) %*% y)

    if (use_center(x_seed)) {
        out <- sweep(out, 2, as.numeric(get_center(x_seed) %*% y), "-", check.margin=FALSE)
    }

    out
}

#' @export
#' @importFrom Matrix t 
#' @importFrom DelayedArray seed DelayedArray
setMethod("%*%", c("ANY", "ScaledMatrix"), function(x, y) {
    y_seed <- seed(y)
    if (is_transposed(y_seed)) {
        if (!is.null(dim(x))) {
            # Vectors don't quite behave as 1-column matrices here.
            # so we need to be a bit more careful.
            x <- t(x) 
        }
        out <- t(.rightmult_ScaledMatrix(y_seed, x))
    } else {
        out <- .leftmult_ScaledMatrix(x, y_seed)
    }
    DelayedArray(out)
})

#' @importFrom Matrix rowSums
#' @importFrom DelayedArray sweep
.leftmult_ScaledMatrix <- function(x, y_seed) { 
    out <- as.matrix(x %*% get_matrix2(y_seed))

    if (use_center(y_seed)) {
        if (is.null(dim(x))) {
            out <- out - get_center(y_seed) * sum(x)
        } else {
            out <- out - outer(rowSums(x), get_center(y_seed), "*")
        }
    }

    if (use_scale(y_seed)) {
        out <- sweep(out, 2, get_scale(y_seed), "/", check.margin=FALSE)
    }

    out
}

#' @export
#' @importFrom DelayedArray seed DelayedArray
setMethod("%*%", c("ScaledMatrix", "ScaledMatrix"), function(x, y) {
    x_seed <- seed(x)
    y_seed <- seed(y)
    res <- .dual_mult_dispatcher(x_seed, y_seed, is_transposed(x_seed), is_transposed(y_seed))
    DelayedArray(res)
})

#' @importFrom Matrix t
.dual_mult_dispatcher <- function(x_seed, y_seed, x_trans, y_trans) {
    if (!x_trans) {
        if (!y_trans) {
            res <- .multiply_u2u(x_seed, y_seed)
        } else {
            res <- .multiply_u2t(x_seed, y_seed)
        }
    } else {
        if (!y_trans) {
            res <- .multiply_t2u(x_seed, y_seed)
        } else {
            res <- .multiply_u2u(y_seed, x_seed)
            res <- t(res)
        }
    }
    res
}

###################################
# ScMat %*% ScMat utilities.

# We do not implement ScMat %*% ScMat in terms of left/right %*%.
# This would cause scaling to be applied on one of the ScMats,
# collapsing it into a DelayedMatrix. Subsequent multiplication 
# would use block processing, which would be too slow.

#' @importFrom Matrix drop rowSums
#' @importFrom DelayedArray sweep
.multiply_u2u <- function(x_seed, y_seed) 
# Considering the problem of (X - C_x)S_x (Y - C_y)S_y.
{
    # Computing X S_x Y S_y
    x0 <- get_matrix2(x_seed)
    if (use_scale(x_seed)) {
        x0 <- ScaledMatrix(x0, scale=get_scale(x_seed))
    } 

    result <- as.matrix(x0 %*% get_matrix2(y_seed))
    if (use_scale(y_seed)) {
        result <- sweep(result, 2, get_scale(y_seed), "/", check.margin=FALSE)
    }

    # Computing C_x S_x Y S_y, and subtracting it from 'result'.
    if (use_center(x_seed)) {
        x.center <- get_center(x_seed)
        if (use_scale(x_seed)) {
            x.center <- x.center / get_scale(x_seed)
        }

        component2 <- drop(x.center %*% get_matrix2(y_seed))
        if (use_scale(y_seed)) {
            component2 <- component2 / get_scale(y_seed)
        }

        result <- sweep(result, 2, component2, "-", check.margin=FALSE)
    }

    # Computing C_x S_x C_y S_y, and adding it to 'result'.
    if (use_center(x_seed) && use_center(y_seed)) {
        x.center <- get_center(x_seed)
        if (use_scale(x_seed)) {
            x.center <- x.center / get_scale(x_seed)
        }

        y.center <- get_center(y_seed)
        if (use_scale(y_seed)) {
            y.center <- y.center / get_scale(y_seed)
        }

        component4 <- sum(x.center) * y.center
        result <- sweep(result, 2, component4, "+", check.margin=FALSE)
    }

    # Computing X S_x C_y S_y, and subtracting it from 'result'.
    # This is done last to avoid subtracting large values.
    if (use_center(y_seed)) {
        y.center <- get_center(y_seed)
        if (use_scale(y_seed)) {
            y.center <- y.center / get_scale(y_seed)
        }

        component3 <- outer(rowSums(x0), y.center)
        result <- result - component3
    }

    result
}

#' @importFrom Matrix tcrossprod drop
#' @importFrom DelayedArray sweep
.multiply_u2t <- function(x_seed, y_seed) 
# Considering the problem of (X - C_x)S_x S_y(Y' - C_y')
{
    # Computing X S_x S_y Y'
    x0 <- get_matrix2(x_seed)
    if (use_scale(x_seed) || use_scale(y_seed)) {
        scaling <- 1
        if (use_scale(x_seed)) {
            scaling <- scaling * get_scale(x_seed)
        }
        if (use_scale(y_seed)) {
            scaling <- scaling * get_scale(y_seed)
        }
        x0 <- ScaledMatrix(x0, scale=scaling)
    }
    result <- as.matrix(tcrossprod(x0, get_matrix2(y_seed)))

    # Computing C_x S_x S_y Y', and subtracting it from 'result'.
    if (use_center(x_seed)) {
        x.center <- get_center(x_seed)
        if (use_scale(x_seed)) {
            x.center <- x.center / get_scale(x_seed)
        }
        if (use_scale(y_seed)) {
            x.center <- x.center / get_scale(y_seed)
        }

        component2 <- drop(tcrossprod(x.center, get_matrix2(y_seed)))
        result <- sweep(result, 2, component2, "-", check.margin=FALSE)
    }

    # Computing C_x S_x S_y C_y', and adding it to 'result'.
    if (use_center(x_seed) && use_center(y_seed)) {
        x.center <- get_center(x_seed)
        if (use_scale(x_seed)) {
            x.center <- x.center / get_scale(x_seed)
        }

        y.center <- get_center(y_seed)
        if (use_scale(y_seed)) {
            y.center <- y.center / get_scale(y_seed)
        }

        component4 <- sum(x.center*y.center)
        result <- result + component4
    }

    # Computing X S_x S_y C_y', and subtracting it from 'result'.
    # This is done last to avoid subtracting large values.
    if (use_center(y_seed)) {
        component3 <- drop(x0 %*% get_center(y_seed))
        result <- result - component3
    }

    result
}

#' @importFrom Matrix crossprod colSums
#' @importFrom DelayedArray sweep
.multiply_t2u <- function(x_seed, y_seed) 
# Considering the problem of S_x(X' - C_x') (Y - C_y)S_y
{
    # C mputing X' Y 
    x0 <- get_matrix2(x_seed)
    y0 <- get_matrix2(y_seed)
    result <- as.matrix(crossprod(x0, y0))

    # Computing C_x' Y, and subtracting it from 'result'.
    if (use_center(x_seed)) {
        x.center <- get_center(x_seed)
        component2 <- outer(x.center, colSums(y0))
        result <- result - component2
    }

    # Computing C_x' C_y, and adding it to 'result'.
    if (use_center(x_seed) && use_center(y_seed)) {
        x.center <- get_center(x_seed)
        y.center <- get_center(y_seed)
        component4 <- outer(x.center, y.center) * nrow(y0)
        result <- result + component4
    }

    # Computing X' C_y, and subtracting it from 'result'.
    # This is done last to avoid subtracting large values.
    if (use_center(y_seed)) {
        component3 <- outer(colSums(x0), get_center(y_seed))
        result <- result - component3
    }

    if (use_scale(x_seed)) {
        result <- result / get_scale(x_seed)
    } 
    if (use_scale(y_seed)) {
        result <- sweep(result, 2, get_scale(y_seed), "/", check.margin=FALSE)
    }

    result
}

###################################
# Cross-product. 

# Technically, we could implement this in terms of '%*%',
# but we use specializations to exploit native crossprod() for '.matrix',
# which is probably more efficient.

#' @export
#' @importFrom Matrix crossprod
#' @importFrom DelayedArray seed DelayedArray
setMethod("crossprod", c("ScaledMatrix", "missing"), function(x, y) {
    x_seed <- seed(x)
    if (is_transposed(x_seed)) {
        # No need to t(), the output is symmetric anyway. 
        out <- .tcp_ScaledMatrix(x_seed)
    } else {
        out <- .cross_ScaledMatrix(x_seed)
    }

    DelayedArray(out)
})

#' @importFrom Matrix crossprod colSums
#' @importFrom DelayedArray sweep
.cross_ScaledMatrix <- function(x_seed) {
    x0 <- get_matrix2(x_seed)
    out <- as.matrix(crossprod(x0))

    if (use_center(x_seed)) {
        centering <- get_center(x_seed)
        colsums <- colSums(x0)

        # Minus, then add, then minus, to mitigate cancellation.
        out <- out - outer(centering, colsums)
        out <- out + outer(centering, centering) * nrow(x0)
        out <- out - outer(colsums, centering)
    }

    if (use_scale(x_seed)) {
        scaling <- get_scale(x_seed)
        out <- sweep(out / scaling, 2, scaling, "/", check.margin=FALSE)
    }
    out
}

#' @export
#' @importFrom Matrix crossprod
#' @importFrom DelayedArray seed DelayedArray
setMethod("crossprod", c("ScaledMatrix", "ANY"), function(x, y) {
    x_seed <- seed(x)
    if (is_transposed(x_seed)) {
        out <- .rightmult_ScaledMatrix(x_seed, y)
    } else {
        out <- .rightcross_ScaledMatrix(x_seed, y)
    }
    DelayedArray(out)
})

#' @importFrom Matrix crossprod colSums
.rightcross_ScaledMatrix <- function(x_seed, y) {
    out <- as.matrix(crossprod(get_matrix2(x_seed), y))

    if (use_center(x_seed)) {
        if (is.null(dim(y))) {
            out <- out - get_center(x_seed) * sum(y)
        } else {
            out <- out - outer(get_center(x_seed), colSums(y))
        }
    }
    
    if (use_scale(x_seed)) {
        out <- out / get_scale(x_seed)
    }

    out
}

#' @export
#' @importFrom Matrix crossprod
#' @importFrom DelayedArray seed DelayedArray
setMethod("crossprod", c("ANY", "ScaledMatrix"), function(x, y) {
    y_seed <- seed(y)
    if (is_transposed(y_seed)) {
        out <- t(.rightmult_ScaledMatrix(y_seed, x))
    } else {
        out <- .leftcross_ScaledMatrix(x, y_seed)
    }
    DelayedArray(out)
})

#' @importFrom Matrix crossprod colSums
#' @importFrom DelayedArray sweep
.leftcross_ScaledMatrix <- function(x, y_seed) {
    out <- as.matrix(crossprod(x, get_matrix2(y_seed)))

    if (use_center(y_seed)) {
        if (is.null(dim(x))) {
            out <- sweep(out, 2, sum(x) * get_center(y_seed), "-", check.margin=FALSE)
        } else {
            out <- out - outer(colSums(x), get_center(y_seed))
        }
    }

    if (use_scale(y_seed)) {
        out <- sweep(out, 2, get_scale(y_seed), "/", check.margin=FALSE)
    }

    out
}

#' @export
#' @importFrom Matrix crossprod
#' @importFrom DelayedArray DelayedArray seed
setMethod("crossprod", c("ScaledMatrix", "ScaledMatrix"), function(x, y) {
    x_seed <- seed(x)
    y_seed <- seed(y)
    res <- .dual_mult_dispatcher(x_seed, y_seed, !is_transposed(x_seed), is_transposed(y_seed))
    DelayedArray(res)
})

###################################
# Transposed cross-product. 

# Technically, we could implement this in terms of '%*%',
# but we use specializations to exploit native tcrossprod() for '.matrix',
# which is probably more efficient.

#' @export
#' @importFrom Matrix tcrossprod
#' @importFrom DelayedArray seed DelayedArray sweep
setMethod("tcrossprod", c("ScaledMatrix", "missing"), function(x, y) {
    x_seed <- seed(x)
    if (is_transposed(x_seed)) {
        out <- .cross_ScaledMatrix(x_seed)
    } else {
        out <- .tcp_ScaledMatrix(x_seed)
    }
    DelayedArray(out)
})

#' @importFrom Matrix tcrossprod
.tcp_ScaledMatrix <- function(x_seed) {
    x0 <- get_matrix2(x_seed)

    if (use_scale(x_seed)) {
        out <- as.matrix(.internal_tcrossprod(x0, get_scale(x_seed)))
    } else {
        out <- as.matrix(tcrossprod(x0))
    }
    
    if (use_center(x_seed)) {
        centering <- get_center(x_seed)

        if (use_scale(x_seed)) {
            centering <- centering / get_scale(x_seed)
            extra <- centering / get_scale(x_seed)
        } else {
            extra <- centering
        }
            
        # With scaling, the use of 'extra' mimics sweep(x0, 2, get_scale(x), "/"),
        # except that the scaling is applied to 'centering' rather than directly to 'x0'.
        # Without scaling, 'extra' and 'centering' are interchangeable.
        component <- tcrossprod(extra, x0)

        # Minus, then add, then minus, to mitigate cancellation.
        out <- sweep(out, 2, as.numeric(component), "-", check.margin=FALSE)
        out <- out + sum(centering^2)
        out <- out - as.numeric(x0 %*% extra)
    }

    out
}

#' @export
#' @importFrom Matrix tcrossprod t
#' @importFrom DelayedArray seed DelayedArray sweep
setMethod("tcrossprod", c("ScaledMatrix", "ANY"), function(x, y) {
    if (is.null(dim(y))) { # for consistency with base::tcrossprod.
        stop("non-conformable arguments")
    }

    x_seed <- seed(x)
    if (is_transposed(x_seed)) {
        out <- t(.leftmult_ScaledMatrix(y, x_seed))
    } else {
        out <- .righttcp_ScaledMatrix(x_seed, y)
    }
    DelayedArray(out)
})

#' @importFrom Matrix tcrossprod
.righttcp_ScaledMatrix <- function(x_seed, y) {
    if (use_scale(x_seed)) {
        # 'y' cannot be a vector anymore, due to the check above.
        y <- sweep(y, 2, get_scale(x_seed), "/", check.margin=FALSE)
    }

    out <- as.matrix(tcrossprod(get_matrix2(x_seed), y))

    if (use_center(x_seed)) {
        out <- sweep(out, 2, as.numeric(tcrossprod(get_center(x_seed), y)), "-", check.margin=FALSE)
    }

    out
}

#' @export
#' @importFrom Matrix tcrossprod t
#' @importFrom DelayedArray seed DelayedArray
setMethod("tcrossprod", c("ANY", "ScaledMatrix"), function(x, y) {
    y_seed <- seed(y) 
    if (is_transposed(y_seed)) {
        out <- .leftmult_ScaledMatrix(x, y_seed)
    } else {
        out <- .lefttcp_ScaledMatrix(x, y_seed)
    }
    DelayedArray(out)
})

#' @importFrom Matrix tcrossprod 
.lefttcp_ScaledMatrix <- function(x, y_seed) {
    if (use_scale(y_seed)) {
        if (is.null(dim(x))) {
            x <- x / get_scale(y_seed)
        } else { 
            x <- sweep(x, 2, get_scale(y_seed), "/", check.margin=FALSE)
        }
    }

    out <- as.matrix(tcrossprod(x, get_matrix2(y_seed)))

    if (use_center(y_seed)) {
        out <- out - as.numeric(x %*% get_center(y_seed))
    }
    out
}

#' @export
#' @importFrom Matrix tcrossprod
#' @importFrom DelayedArray DelayedArray seed
setMethod("tcrossprod", c("ScaledMatrix", "ScaledMatrix"), function(x, y) {
    x_seed <- seed(x)
    y_seed <- seed(y)
    res <- .dual_mult_dispatcher(x_seed, y_seed, is_transposed(x_seed), !is_transposed(y_seed))
    DelayedArray(res)
})

###################################
# Extra code for corner-case calculations of the transposed cross-product.

#' @importFrom DelayedArray seed DelayedArray
.update_scale <- function(x, s) {
    x_seed <- seed(x)
    if (use_scale(x_seed)) {
        s <- s * get_scale(x_seed)
    }
    x_seed@scale <- s
    x_seed@use_scale <- TRUE
    DelayedArray(x_seed)
}

#' @importFrom Matrix tcrossprod 
#' @importFrom methods is
#' @importFrom DelayedArray seed
.internal_tcrossprod <- function(x, scale.) 
# Computes tcrossprod(sweep(x, 2, scale, "/")) when 'x' is a matrix-like object.
# 'scale' can be assumed to be non-NULL here.
# This will always return a dense ordinary matrix.
{
    if (!is(x, "ScaledMatrix")) {
        x <- sweep(x, 2, scale., "/", check.margin=FALSE) 
        return(as.matrix(tcrossprod(x)))
    }

    x_seed <- seed(x)
    if (!is_transposed(x_seed)) {
        x <- .update_scale(x, scale.)
        return(as.matrix(tcrossprod(x)))
    }

    inner <- get_matrix2(x_seed)
    if (is(inner, "ScaledMatrix")) {
        if (is_transposed(seed(inner))) {
            component1 <- as.matrix(crossprod(.update_scale(inner, scale.)))
        } else {
            component1 <- .internal_tcrossprod(t(inner), scale.) # recurses. 
        }
    } else {
        component1 <- as.matrix(crossprod(inner/scale.))
    }
           
    if (use_center(x_seed)) {
        centering <- get_center(x_seed)
        component2 <- .internal_mult_special(centering, scale., inner)
        component3 <- t(component2)
        component4 <- outer(centering, centering) * sum(1/scale.^2)
        final <- (component1 - component2) + (component4 - component3)
    } else {
        final <- component1
    }

    if (use_scale(x_seed)) {
        x.scale <- get_scale(x_seed)
        final <- final / x.scale
        final <- sweep(final, 2, x.scale, "/", check.margin=FALSE) 
    }

    final 
}

#' @importFrom methods is
#' @importFrom DelayedArray seed
.internal_mult_special <- function(center, scale., Z)
# Computes C^T * S^2 * Z where C is a matrix of 'centers' copied byrow=TRUE;
# S is a diagonal matrix filled with '1/scale'; and 'Z' is a ScaledMatrix.
# This will always return a dense ordinary matrix.
{
    if (!is(Z, "ScaledMatrix")) {
        return(outer(center, colSums(Z / scale.^2)))
    }
    
    Z_seed <- seed(Z)
    if (is_transposed(Z_seed)) {
        Z <- .update_scale(Z, scale.^2)
        return(outer(center, colSums(Z)))
    }

    output <- .internal_mult_special(center, scale., get_matrix2(Z_seed)) # recurses.

    if (use_center(Z_seed)) {
        output <- output - outer(center, get_center(Z_seed)) * sum(1/scale.^2)
    }

    if (use_scale(Z_seed)) {
        output <- sweep(output, 2, get_scale(Z_seed), "/")
    }

    output
}
