caensvm <- function(x, y, lam1, tau, theta, s, kernel = "linear", ...){
  ## It is implemented based on the Pegasos-based DC algorithm
  ## Pegasos was proposed by Shalev-Shwartz et al. (2011)
  ## all vectors are in column-form
  sgaEN <- function(x, y, v, tau, theta){# sub-gradient of aEN loss function
    xn <- nrow(x)
    xp <- ncol(x)
    sg <- matrix(0, nrow = xp, ncol = 1)
    for (i in 1:xn) {
      cond <- 1 - y[i]*c(x[i,]%*%v)
      if(!is.na(cond)){
        if(cond < 0){
          sg <- sg + tau*(theta*c(x[i,]%*%v)*x[i,] + (1 - 2*theta)*y[i]*x[i,])
        }else{
          sg <- sg + theta*c(x[i,]%*%v)*x[i,] - y[i]*x[i,]
        }
      }
    }
    return(as.matrix(sg))
  }
  sgAEN1 <- function(x, y, w, tau, theta, s){# sub-gradient of aEN1 loss function
    xn <- nrow(x)
    xp <- ncol(x)
    u1 <- (1 - theta - sqrt((1 - theta)^2 + 2*theta*s/tau))/theta
    u2 <- (theta - 1 + sqrt((1 - theta)^2 + 2*theta*s))/theta
    sg <- matrix(0, nrow = xp, ncol = 1)
    for (i in 1:xn) {
      cond <- 1 - y[i]*c(x[i,]%*%w)
      if(cond < u1){
        sg <- sg + tau*(theta*c(x[i,]%*%w)*x[i,] + (1 - 2*theta)*y[i]*x[i,])
      }else if(cond > u2){
        sg <- sg + theta*c(x[i,]%*%w)*x[i,] - y[i]*x[i,]
      }
    }
    return(as.matrix(sg))
  }
  pegasos <- function(x, y, w, m, lam1, tau, theta, s){
    t2 <- 0
    v <- w
    while(t2 <= 500){
      t2 <- t2 + 1
      At <- sample(nrow(x), m)
      xm <- x[At,]
      ym <- y[At]
      
      # main iteration part if Pegasos
      dF <- v + (lam1/m)*(sgaEN(xm, ym, v, tau, theta) - sgAEN1(xm, ym, w, tau, theta, s))
      v <- v - (lam1/t2)*dF
    }
    return(v)
  }
  
  # main part of Pegasos-based DC algorithm for CaENSVM
  xn <- nrow(x)
  if(kernel == "gaussian"){
    randx <- x[sample(xn, floor(0.2*xn)),]
    x <- gaussian_kernel(x, t(randx), ...)
  }
  xp <- ncol(x)
  xe <- matrix(1, nrow = xn, ncol = 1)
  xx <- cbind(x, xe)
  w0 <- matrix(0, nrow = xp+1, ncol = 1)
  w1 <- pegasos(x = xx, y, w = w0, m = 10, lam1, tau, theta, s)
  w1[is.na(w1)] <- 1
  
  t1 <- 0
  while(t1 < 10 && sqrt(sum((w0 - w1)^2)) >= 1e-3){
    t1 <- t1 + 1
    w0 <- w1
    w1 <- pegasos(x = xx, y, w = w0, m = 10, lam1, tau, theta, s)
    w1[is.na(w1)] <- 1
  }
  
  wnorm <- sqrt(sum(w1[1:xp]^2))
  if(kernel == "linear"){
    caenClassifier <- list(w = as.matrix(w1[1:xp]/wnorm), b = w1[xp+1]/wnorm)
  }else if(kernel == "gaussian"){
    caenClassifier <- list(w = as.matrix(w1[1:xp]/wnorm), b = w1[xp+1]/wnorm, randx = randx)
  }
  return(caenClassifier)
}
