#library("quantreg")
#library("MASS")
#library("pracma") #gramSchmidt

#HPCA(Method="E");IQR()
Fregularize <- function(F){
  Freg <- Orthogonalize(F)
  return(Freg)
}

Lregularize <- function(L){
  svdL  <- svd(L)
  svdLD <- diag(svdL$d)
  Lreg <- (svdL$u)%*%svdLD
  return(Lreg)
}

Orthogonalize <- function(Z){
  gs <- gramSchmidt(Z)
  Q <- gs$Q
  return(Q)
}
