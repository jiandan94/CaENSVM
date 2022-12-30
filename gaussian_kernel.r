gaussian_kernel <- function(A, B, sigma) {
    A <- as.matrix(A)
    B <- as.matrix(B)
    n <- nrow(A)
    p <- ncol(B)

    repeatMatrix <- function(x, n, p, ...){
        # setdiff(x, y): find the elements in x but not in y
        if(n == p){
            repnum <- n
        }else{
            repnum <- setdiff(c(n, p), length(x))
        }
        return(matrix(rep(x, repnum), ...))
    }

    if(ncol(A) == nrow(B)){
        rowAnorm <- apply(A, 1, function(x){return(sum(x^2))})
        AnormMatrix <- repeatMatrix(rowAnorm, n, p, nrow = n, byrow = F)

        colBnorm <- apply(B, 2, function(x){return(sum(x^2))})
        BnormMatrix <- repeatMatrix(colBnorm, n, p, ncol = p, byrow = T)

        K <- exp(-sigma*(AnormMatrix + BnormMatrix - 2*A%*%B))
    }else{
        return("ERROR: non-conformable arguments")
    }
    return(K)
}