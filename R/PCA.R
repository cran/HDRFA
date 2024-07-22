PCA=function(X,r,constraint="L"){
  p=ncol(X);n=nrow(X)

  X_svd=svd(X,r,r)

  if(constraint=="L"){
    Lhat=sqrt(p)*X_svd$v[,1:r]
    Fhat=X%*%Lhat/p
    return(list(Fhat=as.matrix(Fhat),Lhat=as.matrix(Lhat)))
  }

  if(constraint=="F"){
    Fhat=sqrt(n)*X_svd$u[,1:r]
    Lhat=t(X)%*%Fhat/n
    return(list(Fhat=as.matrix(Fhat),Lhat=as.matrix(Lhat)))
  }
}

PCA_FN=function(X,rmax){
  rmax=rmax+1
  p=ncol(X);n=nrow(X)

  X_svd=svd(X,rmax,rmax)
  rhat=which.max(X_svd$d[1:(rmax-1)]/(X_svd$d[2:rmax]))
  return(rhat)
}
