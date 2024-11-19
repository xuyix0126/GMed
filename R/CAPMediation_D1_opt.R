#################################################
CAPMediation_D1_opt<-function(X,M,Y,H=NULL,max.itr=1000,tol=1e-4,trace=FALSE,score.return=TRUE,theta0.mat=NULL,ninitial=NULL,seed=100)
{
  # X: (X, W) matrix
  # M: a list of length n, M
  # Y: n by 1 vector
  # H: a p by p positive definite matrix for theta constraint

  n<-length(Y)
  p<-ncol(M[[1]])

  nT<-rep(NA,n)
  M.cov<-array(NA,c(p,p,n))
  for(i in 1:n)
  {
    nT[i]<-nrow(M[[i]])
    M.cov[,,i]<-cov(M[[i]])
  }

  re<-CAPMediation_D1_opt_Mcov(X=X,M.cov=M.cov,Y=Y,nT=nT,H=H,max.itr=max.itr,tol=tol,trace=trace,score.return=score.return,theta0.mat=theta0.mat,ninitial=ninitial,seed=seed)

  return(re)
}
#################################################
