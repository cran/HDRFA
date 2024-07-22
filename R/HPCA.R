SC=function(mu_update,mu_init,tol=0.001){

  mu_inl2=norm(mu_init,"F")
  if(mu_inl2==0){
    res=norm((mu_update-mu_init),"F")
  } else{
    res=norm((mu_update-mu_init),"F")/mu_inl2
  }
  return(res<=tol)
}

HPCA=function(X,r,Method="E",tau=NULL,scale_est="MAD",L_init=NULL,F_init=NULL,maxiter_HPCA=100,maxiter_HLM=100,eps=0.001){

  dimen=dim(X)
  T=dimen[1];N=dimen[2]

  if(Method=="E"){
    if(is.null(tau)){
      tau=1.345
    }
    
    if(is.null(L_init) & is.null(F_init)){
      
      fit_init=PCA(X,r)
      F_init=fit_init$Fhat;L_init=fit_init$Lhat
    }

    if(is.null(L_init)){
      L_init=matrix(0,N,r)
      for(i in 1:N){
        L_init[i,]=myrlm(F_init,X[,i],tau,scale_est=scale_est,maxiter=maxiter_HLM)$coefficients
        #L_init[i,]=rlm(F_init,X[,i],maxit=maxiter_HLM)$coefficients
      }
      FLR_fit=FLregularize(L_init,F_init,N,T,r);F_init=FLR_fit$Fhat;L_init=FLR_fit$Lhat
    }

    if(is.null(F_init)){
      F_init=matrix(0,T,r)
      for(t in 1:T){
        F_init[t,]=myrlm(L_init,X[t,],tau,scale_est=scale_est,maxiter=maxiter_HLM)$coefficients
        #F_init[t,]=rlm(L_init,X[t,],maxit=maxiter_HLM)$coefficients
      }
      FLR_fit=FLregularize(L_init,F_init,N,T,r);F_init=FLR_fit$Fhat;L_init=FLR_fit$Lhat
    }
    CC_init=F_init%*%t(L_init)

    F_update=F_init;L_update=L_init;m=0
    while(TRUE){

      for(t in 1:T){
        F_update[t,]=myrlm(L_init,X[t,],tau,scale_est=scale_est,maxiter=maxiter_HLM)$coefficients
        #F_update[t,]=rlm(L_init,X[t,],maxit=maxiter_HLM)$coefficients
      }

      for(i in 1:N){
        L_update[i,]=myrlm(F_update,X[,i],tau,scale_est=scale_est,maxiter=maxiter_HLM)$coefficients
        #L_update[i,]=rlm(F_update,X[,i],maxit=maxiter_HLM)$coefficients
      }

      CC_update=F_update%*%t(L_update)

      if(m>=maxiter_HPCA){
        warning(gettextf("'HPCA' failed to converge in %d steps", maxiter_HPCA))
        
        FLR_fit=FLregularize(L_update,F_update,N,T,r);F_update=FLR_fit$Fhat;L_update=FLR_fit$Lhat
        return(list(Fhat=as.matrix(F_update),Lhat=as.matrix(L_update),iter=m))
      }

      if(SC(CC_update,CC_init,eps)){
        FLR_fit=FLregularize(L_update,F_update,N,T,r);F_update=FLR_fit$Fhat;L_update=FLR_fit$Lhat
        return(list(Fhat=as.matrix(F_update),Lhat=as.matrix(L_update),iter=m))
      } else{
        L_init=L_update
        F_init=F_update
        CC_init=CC_update
        #tau=NULL
      }
      m=m+1
    }
  }

  if(Method=="P"){
    if(is.null(L_init) & is.null(F_init)){
      fit_init=PCA(X,r)
      L_init=fit_init$Lhat
      F_init=fit_init$Fhat
    }

    if(is.null(L_init)){
      fit_init=PCA(X,r)
      L_init=fit_init$Lhat
    }

    if(is.null(F_init)){
      fit_init=PCA(X,r)
      F_init=fit_init$Fhat
    }

    Ehat_init=X-F_init%*%t(L_init)

    m=1
    Loss=c()

    while(TRUE){
      Ehat_l2=apply(Ehat_init,1,norm,type="2")

      if(is.null(tau)){
        Ehat_order=order(Ehat_l2);omega=c();T1=floor(T/2);tau=Ehat_l2[Ehat_order[T1]]

        omega[Ehat_order[1:T1]]=0.5
        for(t in (T1+1):T){
          frac2=sum((X[Ehat_order[t],])^2)-sum((t(X[Ehat_order[t],])%*%L_init)^2)/N
          omega[Ehat_order[t]]=tau/2/sqrt(frac2)
        }

        Sigmahat=matrix(0,N,N)
        for(t in 1:T){
          Sigmahat=Sigmahat+omega[t]*X[t,]%*%t(X[t,])
        }
        Sigmahat=Sigmahat/T
      }

      L_update=sqrt(N)*eigen(Sigmahat)$vectors[,1:r]
      F_update=X%*%L_update/N
      Ehat_update=X-F_update%*%t(L_update)

      if(m>=maxiter_HPCA){
        warning(gettextf("'HPCA' failed to converge in %d steps", maxiter_HPCA))
        return(list(Fhat=as.matrix(F_update),Lhat=as.matrix(L_update),iter=m))
      }

      if(SC(Ehat_update,Ehat_init,eps)){
        return(list(Fhat=as.matrix(F_update),Lhat=as.matrix(L_update),iter=m))
      } else{
        L_init=L_update
        F_init=F_update
        Ehat_init=Ehat_update
        tau=NULL
      }
      m=m+1
    }
  }
}

HPCA_FN=function(X,rmax,Method="E",tau=NULL,scale_est="MAD",threshold=NULL,L_init=NULL,F_init=NULL,maxiter_HPCA=100,maxiter_HLM=100,eps=0.001){

  dimen=dim(X)
  T=dimen[1];N=dimen[2]

  fit=HPCA(X=X,r=rmax,Method=Method,tau=tau,scale_est=scale_est,L_init=L_init,F_init=F_init,maxiter_HPCA=maxiter_HPCA,maxiter_HLM=maxiter_HLM,eps=eps)
  VK=sort(diag((t(fit$Fhat)%*%fit$Fhat)/T),decreasing=TRUE)
  
  if(is.null(threshold)){
    threshold=VK[1]*(min(N,T))^(-1/3)
  }
  rhat=sum(VK>threshold)
  return(rhat)
}
