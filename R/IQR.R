IQR <- function(X,r,tau,max_iter=100,eps=10^(-6))#X n*p T*N  n=T, p=N
{
  p <- ncol(X)
  n <- nrow(X)
  L0hat <- matrix(0,p,r)#p*r
  F0hat=matrix(rnorm(n*r,0,1),n,r)
  #F0hat <- PCA(X,r,constraint="F")$Fhat #the PCA estimator
  F0hat <- Fregularize(F0hat)
  for (i in 1:p){
    L0hat[i,] <- rq(X[,i]~F0hat-1,tau=tau)$coefficients
  }
  Lthat <- Lregularize(L0hat)
  Fthat <- F0hat

  t <- 0
  CCdiff <- 10^6
  while (((CCdiff>eps)&(t<max_iter))){
    CC0=Lthat%*%t(Fthat)
    for (i in 1:n){
      Fthat[i,] <- rq(X[i,]~Lthat-1,tau=tau)$coefficients
    }
    Fthat <- Fregularize(Fthat)
    for (i in 1:p){
      Lthat[i,] <- rq(X[,i]~Fthat-1,tau=tau)$coefficients
    }
    Lthat <- Lregularize(Lthat)
    CC1=Lthat%*%t(Fthat)
    CCdiff=mean(abs(CC1-CC0)/abs(CC0))
    t <- t+1
    if (t==max_iter){
      warning(gettextf("'IQR' failed to converge in %d steps", max_iter))
    }
  }
  return(list(Fhat=Fthat,Lhat=Lthat,t=t))
}

IQR_FN<-function(X,rmax,tau,threshold=NULL,max_iter=100,eps=10^(-6)){

  N <- ncol(X)
  T <- nrow(X)
  rip <- IQR(X,rmax,tau,max_iter=max_iter,eps=eps)
  Lhat <- rip$Lhat
  VK <- eigen((t(Lhat)%*%Lhat)/(N),only.values=TRUE)$values
  if(is.null(threshold)){
    threshold<- VK[1]*(min(N,T))^(-1/3)
  }
  rhat <- sum(VK>threshold)
  return(rhat)
}
