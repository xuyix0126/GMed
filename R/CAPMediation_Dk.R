#################################################
# higher order components
# option:
# 1. remove identified components from M direct
# 2. remove identified components from M ~ X + W residuals

CAPMediation_Dk<-function(X,M,Y,H=NULL,Theta0=NULL,Y.remove=TRUE,max.itr=1000,tol=1e-4,trace=FALSE,score.return=TRUE,theta0.mat=NULL,ninitial=NULL,seed=100)
{
  # X: (X, W) matrix
  # M: a list of length n, M
  # Y: n by 1 vector
  # H: a p by p positive definite matrix for theta constraint
  # Theta0: p by k matrix, identified components

  n<-length(Y)
  p<-ncol(M[[1]])
  nT<-sapply(M,nrow)

  if(is.null(Theta0))
  {
    re<-CAPMediation_D1_opt(X=X,M=M,Y=Y,H=H,max.itr=max.itr,tol=tol,trace=trace,score.return=score.return,theta0.mat=theta0.mat,ninitial=ninitial,seed=seed)
  }else
  {
    p0<-ncol(Theta0)

    M.cov<-array(NA,c(p,p,n))
    for(i in 1:n)
    {
      M.cov[,,i]<-cov(M[[i]])
    }
    score<-matrix(NA,n,p0)
    colnames(score)<-paste0("C",1:p0)
    for(kk in 1:p0)
    {
      score[,kk]<-apply(M.cov,c(3),function(x){return((t(Theta0[,kk])%*%x%*%Theta0[,kk])[1,1])})
    }

    # estimate alpha0
    alpha0.est<-rep(NA,p0)
    alpha.est<-matrix(NA,ncol(X),p0)
    rownames(alpha.est)<-colnames(X)
    colnames(alpha.est)<-paste0("C",1:p0)
    beta.est<-rep(NA,p0)
    alpha0.rnd.est<-matrix(NA,n,p0)
    colnames(alpha0.rnd.est)<-paste0("C",1:p0)
    for(kk in 1:p0)
    {
      otmp<-CAPMediation_coef(X,M,Y,theta=Theta0[,kk])
      alpha0.est[kk]<-otmp$alpha0
      alpha0.rnd.est[,kk]<-otmp$alpha0.rnd
      alpha.est[,kk]<-otmp$alpha
      beta.est[kk]<-otmp$beta
    }

    Mnew<-vector("list",length=n)
    names(Mnew)<-names(M)
    for(i in 1:n)
    {
      Mtmp<-M[[i]]-M[[i]]%*%Theta0%*%t(Theta0)
      Mtmp.svd<-svd(Mtmp)
      # dnew<-c(Mtmp.svd$d[1:(p-p0)],sqrt(exp(alpha0.rnd.est[i,])*nT[i]))
      dnew<-c(Mtmp.svd$d[1:(p-p0)],sqrt(exp(alpha0.rnd.est[i,]-alpha0.est)*nT[i]))
      Mnew[[i]]<-Mtmp.svd$u%*%diag(dnew)%*%t(Mtmp.svd$v)
    }

    if(Y.remove)
    {
      Ynew<-Y-log(score)%*%beta.est
    }else
    {
      Ynew<-Y
    }

    re<-CAPMediation_D1_opt(X=X,M=Mnew,Y=Ynew,H=H,max.itr=max.itr,tol=tol,trace=trace,score.return=score.return,theta0.mat=theta0.mat,ninitial=ninitial,seed=seed)
    re$orthogonal<-c(t(re$theta)%*%Theta0)
  }

  return(re)
}
#################################################
