rwish=function (v, S) 
{
    if (!is.matrix(S)) 
        S <- matrix(S)
    if (nrow(S) != ncol(S)) {
        stop(message = "S not square in rwish().\n")
    }
    if (v < nrow(S)) {
        stop(message = "v is less than the dimension of S in rwish().\n")
    }
    p <- nrow(S)
    CC <- chol(S)
    Z <- matrix(0, p, p)
    diag(Z) <- sqrt(rchisq(p, v:(v - p + 1)))
    if (p > 1) {
        pseq <- 1:(p - 1)
        Z[rep(p * pseq, pseq) + unlist(lapply(pseq, seq))] <- rnorm(p * 
            (p - 1)/2)
    }
    return(crossprod(Z %*% CC))
}

ergMean=function (x, m = 1) 
{
    if (hasArg(m)) 
        m <- max(1, round(m))
    n <- NROW(x)
    if (m > n) 
        stop("Need m <= n")
    if (is.null(dm <- dim(x)) || dm[2] == 1) {
        if (m == 1) 
            ans <- cumsum(x)/1:n
        else if (m == n) 
            ans <- mean(x)
        else ans <- cumsum(c(sum(x[1:m]), x[(m + 1):n]))/m:n
    }
    else {
        if (length(dm) == 2) {
            nm <- colnames(x)
            if (is.null(nm)) 
                nm <- paste(deparse(substitute(x)), 1:NCOL(x), 
                  sep = ".")
            if (m == 1) 
                ans <- apply(x, 2, cumsum)/1:n
            else if (m == n) 
                ans <- matrix(colMeans(x), nrow = 1)
            else ans <- apply(rbind(colSums(x[1:m, ]), x[(m + 
                1):n, ]), 2, cumsum)/m:n
            colnames(ans) <- nm
        }
        else stop("'x' must be a vector or a matrix")
    }
    return(ans)
}

diwish=function (W, v, S) 
{
    if (!is.matrix(S)) 
        S <- matrix(S)
    if (nrow(S) != ncol(S)) {
        stop("S not square in diwish().\n")
    }
    if (!is.matrix(W)) 
        W <- matrix(W)
    if (nrow(W) != ncol(W)) {
        stop("W not square in diwish().\n")
    }
    if (nrow(S) != ncol(W)) {
        stop("W and X of different dimensionality in diwish().\n")
    }
    if (v < nrow(S)) {
        stop("v is less than the dimension of S in  diwish().\n")
    }
    p <- nrow(S)
    gammapart <- sum(lgamma((v + 1 - 1:p)/2))
    ldenom <- gammapart + 0.5 * v * p * log(2) + 0.25 * p * (p - 
        1) * log(pi)
    cholS <- chol(S)
    cholW <- chol(W)
    halflogdetS <- sum(log(diag(cholS)))
    halflogdetW <- sum(log(diag(cholW)))
    invW <- chol2inv(cholW)
    exptrace <- sum(S * invW)
    lnum <- v * halflogdetS - (v + p + 1) * halflogdetW - 0.5 * 
        exptrace
    lpdf <- lnum - ldenom
    return(exp(lpdf))
}

eps <- .Machine$double.eps^0.3

is.pos=function (x) 
{
    if (!is.matrix(x)) 
        stop("x is not a matrix.")
    if (!is.square.matrix(x)) 
        stop("x is not a square matrix.")
    if (!is.symmetric.matrix(x)) 
        stop("x is not a symmetric matrix.")
    eigs <- eigen(x, symmetric = TRUE)$values
		if (any(eigs < 0)) {
	ll=length(which(eigs<0))
	  id=which(eigs<=0) 
	  eigs[id] <- eps
                            }	
    if (any(is.complex(eigs))) 
        return(FALSE)
    if (all(eigs > 0)){ 
        pd <- TRUE
		 }   else {
		 pd <- FALSE}
    return(pd)
}

rmatrix=function (M, U, V) 
{
    if (missing(M)) 
        stop("Matrix M is missing.")
    if (!is.matrix(M)) 
        M <- matrix(M)
    if (missing(U)) 
        stop("Matrix U is missing.")
    if (!is.matrix(U)) 
        U <- matrix(U)
    if (!is.pos(U)) 
        stop("Matrix U is not positive-definite.")
    if (missing(V)) 
        stop("Matrix V is missing.")
    if (!is.matrix(V)) 
        V <- matrix(V)
    if (!is.pos(V)) 
        stop("Matrix V is not positive-definite.")
    if (nrow(M) != nrow(U)) 
        stop("Dimensions of M and U are incorrect.")
    if (ncol(M) != ncol(V)) 
        stop("Dimensions of M and V are incorrect.")
    n <- nrow(U)
    k <- ncol(V)
    Z <- matrix(rnorm(n * k), n, k)
    X <- M + t(chol(U)) %*% Z %*% chol(V)
    return(X)
}
