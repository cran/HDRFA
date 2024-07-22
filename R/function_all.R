FLregularize=function(L,F,N,T,r){
  
  L_svd=svd(L,r,r)
  L_Sigma1=L_svd$v%*%diag(L_svd$d^{-1},r)%*%t(L_svd$v)*sqrt(N)
  L_Sigma2=L_svd$v%*%diag(L_svd$d,r)%*%t(L_svd$v)/sqrt(N)
  F_Sigma=t(F)%*%F/T
  GM=L_Sigma2%*%F_Sigma%*%L_Sigma2
  Gamma=svd(GM)$u
  H=L_Sigma1%*%Gamma
  
  # L_Sigma=t(L)%*%L/N;L_ED=eigen(L_Sigma)
  # L_Sigma1=L_ED$vectors%*%diag(L_ED$values^{-0.5},r)%*%t(L_ED$vectors)
  # L_Sigma2=L_ED$vectors%*%diag(L_ED$values^{0.5},r)%*%t(L_ED$vectors)
  # F_Sigma=t(F)%*%F/T
  # GM=L_Sigma2%*%F_Sigma%*%L_Sigma2
  # Gamma=eigen(GM)$vectors
  # H=L_Sigma1%*%Gamma
  
  L=L%*%H
  F=F%*%solve(t(H))
  
  return(list(Fhat=as.matrix(F),Lhat=as.matrix(L)))
}

