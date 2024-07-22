myrlm=function(X,y,tau=1.345,loss="huber",scale_est="MAD",beta_init=NULL,w_init=NULL,maxiter=100,tol=1e-8){
  
  if(is.null(w_init)){
    w=rep(1,length(y))
  }
  
  if(is.null(beta_init)){
    beta_init=lm.wfit(X,y,w)$coefficients
  }
  
  epsilon_init=y-X%*%beta_init
  
  m=0;Loss=c()
  while(TRUE){
    
    sigmahat=Robust_sigma(epsilon_init,scale_est)
    w=Psi_Type(epsilon_init/sigmahat,tau,loss)
    
    beta_update=lm.wfit(X,y,w)$coefficients
    Loss[m+1]=sum((beta_update-beta_init)^2)/max(1e-20,sum(beta_init^2))
    
    if(Loss[m+1]<tol){
      return(list(coefficients=beta_update,m=m,Loss=Loss))
    } else{
      beta_init=beta_update
      epsilon_init=y-X%*%beta_update
      m=m+1
    }
    
    if(m>maxiter){
      warning(gettextf("'rlm' failed to converge in %d steps", maxiter))
      return(list(coefficients=beta_update,m=m,Loss=Loss))
    }
  }
}

Robust_sigma=function(x,type="MAD"){
  if(type=="MAD"){
    sigmahat=median(abs(x))/0.6745
  }
  
  if(type=="const"){
    sigmahat=1
  }
  return(sigmahat)
}

Psi_Type=function(x,tau,loss){
  if(loss=="huber")  return(pmin(1, tau/abs(x)))
}


