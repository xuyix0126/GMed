# given an estimate of the covariance matrices
CAPMediation_coef_Mcov<-function(X,M.cov,Y,theta)
{
  # X: (X, W) matrix
  # M.cov: p by p by n array, covariance matrix of M
  # Y: n by 1 vector
  # theta: projection vector

  n<-length(Y)
  p<-ncol(M.cov[,,1])

  q<-ncol(X)-1
  if(is.null(colnames(X)))
  {
    if(q>0)
    {
      colnames(X)<-c("X",paste0("W",1:q))
    }else
    {
      colnames(X)<-c("X")
    }
  }

  score<-apply(M.cov,c(3),function(x){return((t(theta)%*%x%*%theta)[1,1])})

  # beta and gamma estimate
  Z<-cbind(rep(1,n),X,log(score))
  colnames(Z)<-c("Intercept",colnames(X),"M")
  mu.est<-c(ginv(t(Z)%*%Z)%*%t(Z)%*%Y)

  gamma0.est<-mu.est[1]
  gamma.est<-mu.est[2:(length(mu.est)-1)]
  beta.est<-mu.est[length(mu.est)]

  # sigma2 estimate
  sigma2.est<-mean((Y-Z%*%mu.est)^2,na.rm=TRUE)

  # estimate alpha0, alpha0.rnd, alpha using mixed effects model
  dtmp<-data.frame(ID=1:n,score=log(score),X)
  eval(parse(text=paste0("fit.tmp<-lme(score~",paste(colnames(X),collapse="+"),",random=~1|ID,data=dtmp,control=lmeControl(opt='optim'))")))
  alpha0.est<-fit.tmp$coefficients$fixed[1]
  alpha.est<-fit.tmp$coefficients$fixed[-1]
  alpha0.rnd.est<-c(fit.tmp$coefficients$random$ID+fit.tmp$coefficients$fixed[1])
  tau2.est<-mean((alpha0.rnd.est-alpha0.est)^2,na.rm=TRUE)

  IE.est<-alpha.est[1]*beta.est

  names(alpha0.est)<-NULL
  names(alpha.est)=names(gamma.est)<-colnames(X)
  re<-list(theta=theta,alpha=alpha.est,beta=beta.est,gamma=gamma.est,IE=IE.est,alpha0=alpha0.est,alpha0.rnd=alpha0.rnd.est,gamma0=gamma0.est,tau2=tau2.est,sigma2=sigma2.est)

  return(re)
}

