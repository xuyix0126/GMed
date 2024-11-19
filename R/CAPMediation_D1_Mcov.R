#################################################
# given an estimate of the covariance matrices
CAPMediation_D1_Mcov<-function(X,M.cov,Y,nT,H=NULL,max.itr=1000,tol=1e-4,trace=FALSE,score.return=TRUE,theta0=NULL)
{
  # X: (X, W) matrix
  # M.cov: p by p by n array, covariance matrix of M
  # Y: n by 1 vector
  # nT: n by 1 vector, # of observation of each subject
  # H: a p by p positive definite matrix for theta constraint

  n<-length(Y)

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

  if(is.null(H))
  {
    H<-apply(M.cov,c(1,2),mean,na.rm=TRUE)
  }

  if(is.null(theta0))
  {
    theta0<-rep(1/sqrt(p),p)
  }

  # initial values
  otmp0<-CAPMediation_coef_Mcov(X,M.cov,Y,theta0)
  alpha.ini<-otmp0$alpha
  beta.ini<-otmp0$beta
  gamma.ini<-otmp0$gamma
  alpha0.ini<-otmp0$alpha0
  alpha0.rnd.ini<-otmp0$alpha0.rnd
  gamma0.ini<-otmp0$gamma0
  tau2.ini<-otmp0$tau2
  sigma2.ini<-otmp0$sigma2

  if(trace)
  {
    alpha.trace<-NULL
    beta.trace<-NULL
    gamma.trace<-NULL
    alpha0.trace<-NULL
    alpha0.rnd.trace<-NULL
    gamma0.trace<-NULL
    theta.trace<-NULL
    tau2.trace<-NULL
    sigma2.trace<-NULL

    obj<-NULL
  }

  s<-0
  diff<-100
  while(s<=max.itr&diff>tol)
  {
    s<-s+1

    obj0<-obj.func(X,M.cov,Y,nT,theta0,alpha0.ini,alpha0.rnd.ini,alpha.ini,beta.ini,gamma0.ini,gamma.ini,tau2.ini,sigma2.ini)
    score<-apply(M.cov,c(3),function(x){return((t(theta0)%*%x%*%theta0)[1,1])})

    # update theta
    U<-c(exp(-alpha0.rnd.ini-X%*%alpha.ini))
    V<-c(Y-gamma0.ini-X%*%gamma.ini)
    A<-matrix(0,p,p)
    for(i in 1:n)
    {
      A<-A+(nT[i]*U[i]-(2*beta.ini*(V[i]-beta.ini*log(score[i])))/(sigma2.ini*score[i]))*M.cov[,,i]
    }
    theta.new<-eigen.solve(A,H)

    otmp.new<-CAPMediation_coef_Mcov(X,M.cov,Y,theta.new)
    alpha.new<-otmp.new$alpha
    beta.new<-otmp.new$beta
    gamma.new<-otmp.new$gamma
    alpha0.new<-otmp.new$alpha0
    alpha0.rnd.new<-otmp.new$alpha0.rnd
    gamma0.new<-otmp.new$gamma0
    tau2.new<-otmp.new$tau2
    sigma2.new<-otmp.new$sigma2

    obj.new<-obj.func(X,M.cov,Y,nT,theta.new,alpha0.new,alpha0.rnd.new,alpha.new,beta.new,gamma0.new,gamma.new,tau2.new,sigma2.new)

    diff<-max(abs(c(alpha.new-alpha.ini,beta.new-beta.ini,gamma.new-gamma.ini)))

    alpha.ini<-alpha.new
    beta.ini<-beta.new
    gamma.ini<-gamma.new
    alpha0.ini<-alpha0.new
    alpha0.rnd.ini<-alpha0.rnd.new
    gamma0.ini<-gamma0.new
    tau2.ini<-tau2.new
    sigma2.ini<-sigma2.new
    theta0<-theta.new

    if(trace)
    {
      alpha.trace<-cbind(alpha.trace,alpha.ini)
      beta.trace<-c(beta.trace,beta.ini)
      gamma.trace<-cbind(gamma.trace,gamma.ini)
      alpha0.trace<-cbind(alpha0.trace,alpha0.ini)
      alpha0.rnd.trace<-cbind(alpha0.rnd.trace,alpha0.rnd.ini)
      gamma0.trace<-cbind(gamma0.trace,gamma0.ini)
      tau2.trace<-c(tau2.trace,tau2.ini)
      sigma2.trace<-c(sigma2.trace,sigma2.ini)
      theta.trace<-cbind(theta.trace,theta0)

      obj<-c(obj,obj.new)
    }

    # print(c(diff,obj.new))
  }

  theta.est<-theta0/sqrt(sum(theta0^2))
  if(theta.est[which.max(abs(theta.est))]<0)
  {
    theta.est<--theta.est
  }
  otmp.est<-CAPMediation_coef_Mcov(X,M.cov,Y,theta.est)
  score<-apply(M.cov,c(3),function(x){return((t(theta.est)%*%x%*%theta.est)[1,1])})
  obj.est<-obj.func(X,M.cov,Y,nT,theta.est,alpha0=otmp.est$alpha0,alpha0.rnd=otmp.est$alpha0.rnd,alpha=otmp.est$alpha,beta=otmp.est$beta,gamma0=otmp.est$gamma0,gamma=otmp.est$gamma,
                    tau2=otmp.est$tau2,sigma2=otmp.est$sigma2)

  re<-otmp.est
  if(score.return)
  {
    re$score<-score
  }
  re$obj<-obj.est

  if(trace)
  {
    rownames(alpha.trace)=rownames(gamma.trace)<-colnames(X)

    re$theta.trace<-theta.trace
    re$alpha.trace<-alpha.trace
    re$beta.trace<-beta.trace
    re$gamma.trace<-gamma.trace
    re$alpha0.trace<-alpha0.trace
    re$alpha0.rnd.trace<-alpha0.rnd.trace
    re$gamma0.trace<-gamma0.trace
    re$tau2.trace<-tau2.trace
    re$sigma2.trace<-sigma2.trace
    re$obj.trace<-obj
  }

  return(re)
}
