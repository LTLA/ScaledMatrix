scale_and_center <- function(y, ref, code) {
    center <- scale <- NULL

    if (code==1L) {
        center <- colMeans(ref)
        scale <- runif(ncol(ref))
        ref <- scale(ref, center=center, scale=scale)
    } else if (code==2L) {
        center <- rnorm(ncol(ref))
        ref <- scale(ref, center=center, scale=FALSE)
    } else if (code==3L) {
        scale <- runif(ncol(ref))
        ref <- scale(ref, center=FALSE, scale=scale)
    }

    # Getting rid of excess attributes.
    attr(ref, "scaled:center") <- NULL
    attr(ref, "scaled:scale") <- NULL

    def <- ScaledMatrix(y, center=center, scale=scale)
    list(def=def, ref=ref)
}

spawn_scenarios_basic <- function(NR, NC, CREATOR, REALIZER) {
    output <- vector("list", 8)
    counter <- 1L

    for (trans in c(FALSE, TRUE)) {
        for (it in 1:4) {
            if (trans) { 
                # Ensure output matrix has NR rows and NC columns after t().
                y <- CREATOR(NC, NR)
            } else {
                y <- CREATOR(NR, NC)
            }
            ref <- REALIZER(y) 

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

spawn_scenarios <- function(NR=50, NC=20) {
    c(
        spawn_scenarios_basic(NR, NC,
            CREATOR=function(r, c) {
                matrix(rnorm(r*c), ncol=c)
            }, 
            REALIZER=identity
        ),
        spawn_scenarios_basic(NR, NC, 
            CREATOR=function(r, c) {
                Matrix::rsparsematrix(r, c, 0.1)
            },
            REALIZER=as.matrix
        )
    )
}

expect_equal_product <- function(x, y) {
    expect_s4_class(x, "DelayedMatrix")
    X <- as.matrix(x)

    # standardize NULL dimnames.
    if (all(lengths(dimnames(X))==0L)) dimnames(X) <- NULL 
    if (all(lengths(dimnames(y))==0L)) dimnames(y) <- NULL 
    expect_equal(X, y)
}

purgenames <- function(mat) {
    if (identical(dimnames(mat), list(NULL, NULL))) {
        dimnames(mat) <- NULL
    }
    mat
}
