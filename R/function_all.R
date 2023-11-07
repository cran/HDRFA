#library("quantreg")
#library("MASS")
#library("pracma") #gramSchmidt

#HPCA(Method="E");IQR()
Fregularize <- function(F,r){
  
  if(r>1){
    Freg <- Orthogonalize(F)
  } else{
    Freg <- F/sqrt(sum(F^2))
  }
  return(Freg)
}

Lregularize <- function(L,r){
  
  if(r>1){
    svdL  <- svd(L)
    svdLD <- diag(svdL$d)
    Lreg <- (svdL$u)%*%svdLD
    return(Lreg)
  } else{
    return(L)
  }
}

Orthogonalize <- function(Z){
  gs <- gramSchmidt(Z)
  Q <- gs$Q
  return(Q)
}
