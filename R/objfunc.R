#################################################
# objective function
obj.func<-function(X,Sigma,Y,nT,theta,alpha0,alpha0.rnd,alpha,beta,gamma0,gamma,tau2,sigma2)
{
  # X: (X, W) matrix
  # Sigma: covariance estimate of M
  # Y: n by 1 outcome vector
  # nT: n by 1 vector, # of observation of each subject
  # theta: projection vector
  # alpha0: fixed intercept of M model
  # alpha0.rnd: n by 1 vector, random effect of alpha0
  # alpha: coefficient of M model
  # beta: coefficient of M in Y model
  # gamma0: intercept of Y model
  # gamma: coefficient of Y model
  # tau2: M model error standard deviation
  # sigma2: Y model error standard deviation

  n<-length(nT)

  # t(theta)%*%Sigma%*%theta
  score<-apply(Sigma,c(3),function(x){return((t(theta)%*%x%*%theta)[1,1])})

  ll1<-sum(((alpha0.rnd+X%*%alpha)+score*exp(-alpha0.rnd-X%*%alpha))*nT)/2
  ll2<-sum((Y-gamma0-X%*%gamma-beta*log(score))^2/sigma2+log(sigma2))/2
  ll3<-sum((alpha0.rnd-alpha0)^2/tau2+log(tau2))/2

  ll<-ll1+ll2+ll3
  return(ll)
}
#################################################
