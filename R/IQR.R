IQR <- function(X,r,tau,L_init=NULL,F_init=NULL,max_iter=100,eps=0.001)
{
  dimen=dim(X)
  T=dimen[1];N=dimen[2]
  
  if(is.null(L_init) & is.null(F_init)){
    
    F_init=PCA(t(X),r)$Lhat
    L_init=matrix(0,N,r)
    for (i in 1:N){
      L_init[i,] <- rq(X[,i]~F_init-1,tau=tau)$coefficients
    }
    L_init <- Lregularize(L_init,r)
  }
  
  if(is.null(L_init)){
  
    L_init=matrix(0,N,r)
    for (i in 1:N){
      L_init[i,] <- rq(X[,i]~F_init-1,tau=tau)$coefficients
    }
    L_init=Lregularize(L_init,r)
  }
  
  if(is.null(F_init)){
    F_init=matrix(0,T,r)
    for (t in 1:T){
      F_init[t,] <- rq(X[t,]~L_init-1,tau=tau)$coefficients
    }
    F_init=Fregularize(F_init,r)
  }
  
  CC_init=F_init%*%t(L_init)
  F_update=as.matrix(F_init);L_update=as.matrix(L_init);m=0
  
  while(TRUE){
    
    for (t in 1:T){
      F_update[t,]=rq(X[t,]~L_init-1,tau=tau)$coefficients
    }
    F_update=Fregularize(F_update,r)
    
    for (i in 1:N){
      L_update[i,]=rq(X[,i]~F_update-1,tau=tau)$coefficients
    }
    L_update=Lregularize(L_update,r)
    
    CC_update=F_update%*%t(L_update)
    
    if(m>=max_iter){
      warning(gettextf("'IQR' failed to converge in %d steps", max_iter))
      F_update=F_update*sqrt(T);L_update=L_update/sqrt(T)
      return(list(Fhat=F_update,Lhat=L_update,iter=m))
    }
    
    if(SC(CC_update,CC_init,eps)){
      F_update=F_update*sqrt(T);L_update=L_update/sqrt(T)
      return(list(Fhat=F_update,Lhat=L_update,iter=m))
    } else{
      L_init=L_update
      F_init=F_update
      CC_init=CC_update
    }
    m=m+1
  }
}

IQR_FN<-function(X,rmax,tau,threshold=NULL,L_init=NULL,F_init=NULL,max_iter=100,eps=10^(-6)){

  N <- ncol(X)
  T <- nrow(X)
  rip <- IQR(X,rmax,tau,L_init=L_init,F_init=F_init,max_iter=max_iter,eps=eps)
  Lhat <- rip$Lhat
  VK <- sort(diag((t(Lhat )%*%Lhat )/N),decreasing=TRUE)
  if(is.null(threshold)){
    threshold<- VK[1]*(min(N,T))^(-1/3)
  }
  rhat <- sum(VK>threshold)
  return(rhat)
}
