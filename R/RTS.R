SK <- function(X){
  p <- ncol(X)
  n <- nrow(X)
  TK <- matrix(0,p,p)
  for(i in 1:(n-2)){
    TT <- matrix(rep(X[i,],n-i),n-i,p,byrow = TRUE)-X[(i+1):n,]
    TT <- t(diag(1/diag(TT%*%t(TT)))%*%TT)%*%TT
    TK <- TK+TT
  }
  TT <- X[n-1,]-X[n,]
  TK <- TK+TT%*%t(TT)/sum(TT^2)
  TK <- 2/(n*(n-1))*TK
  return(TK)
}

RTS <- function(X,r){
  p <- ncol(X)
  n <- nrow(X)
  Khat <- SK(X)
  Lhat <- sqrt(p)*as.matrix(eigen(Khat)$vectors[,1:r]) #p*r
  Fhat <- matrix(0,n,r)#n*r
  for (i in 1:n){
    Fhat[i,] <- lm(X[i,]~Lhat-1)$coefficients
  }
  return(list(Fhat=Fhat,Lhat=Lhat))
}

RTS_FN <- function(X,rmax){
  p <- ncol(X);n <- nrow(X)
  Khat <- SK(X)
  Khat_EV=eigen(Khat,only.values=TRUE)
  rhat=which.max(Khat_EV$values[1:rmax]/(Khat_EV$values[2:(rmax+1)]))
  return(rhat)
}