# multiple initial values
CAPMediation_D1_opt_Mcov<-function(X,M.cov,Y,nT,H=NULL,max.itr=1000,tol=1e-4,trace=FALSE,score.return=TRUE,theta0.mat=NULL,ninitial=NULL,seed=100)
{
  # X: (X, W) matrix
  # M.cov: p by p by n array, covariance matrix of M
  # Y: n by 1 vector
  # nT: n by 1 vector, # of observation of each subject
  # H: a p by p positive definite matrix for theta constraint

  n<-nrow(X)

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

  p<-dim(M.cov)[1]

  #----------------------------------
  if(is.null(theta0.mat)==FALSE)
  {
    if(is.null(ninitial))
    {
      ninitial<-min(c(ncol(theta0.mat),10))
    }else
    {
      ninitial<-min(c(ncol(theta0.mat),ninitial))
    }
  }else
    if(is.null(ninitial))
    {
      ninitial<-min(p,10)
    }

  # theta initial
  if(is.null(theta0.mat))
  {
    set.seed(seed)
    theta.tmp<-matrix(rnorm((max(p,ninitial)+1+5)*p,mean=0,sd=1),nrow=p)
    theta0.mat<-apply(theta.tmp,2,function(x){return(x/sqrt(sum(x^2)))})
  }
  set.seed(seed)
  theta0.mat<-matrix(theta0.mat[,sort(sample(1:ncol(theta0.mat),ninitial,replace=FALSE))],ncol=ninitial)
  #----------------------------------

  #----------------------------------
  # try different initial values with the lowest objective function
  re.tmp<-vector("list",ninitial)
  obj<-rep(NA,ninitial)
  for(kk in 1:ninitial)
  {
    try(re.tmp[[kk]]<-CAPMediation_D1_Mcov(X,M.cov,Y,nT=nT,H=H,max.itr=max.itr,tol=tol,trace=trace,score.return=score.return,theta0=theta0.mat[,kk]))

    if(is.null(re.tmp[[kk]])==FALSE)
    {
      theta.unscale<-re.tmp[[kk]]$theta/sqrt(t(re.tmp[[kk]]$theta)%*%H%*%re.tmp[[kk]]$theta)[1,1]

      try(coef.tmp<-CAPMediation_coef_Mcov(X,M.cov,Y,theta=theta.unscale))
      try(obj[kk]<-obj.func(X,M.cov,Y,nT=nT,theta=theta.unscale,alpha0=coef.tmp$alpha0,alpha0.rnd=coef.tmp$alpha0.rnd,alpha=coef.tmp$alpha,beta=coef.tmp$beta,gamma0=coef.tmp$gamma0,gamma=coef.tmp$gamma,
                            tau2=coef.tmp$tau2,sigma2=coef.tmp$sigma2))
    }
  }
  opt.idx<-which.min(obj)
  re<-re.tmp[[opt.idx]]
  #----------------------------------

  return(re)
}


