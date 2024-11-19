# sample covariance of M
CAPMediation_coef<-function(X,M,Y,theta)
{
  # X: (1, X, W) matrix
  # M: a list of length n, M
  # Y: n by 1 vector
  # theta: projection vector

  M.cov<-array(NA,c(p,p,n))
  for(i in 1:n)
  {
    M.cov[,,i]<-cov(M[[i]])
  }

  re<-CAPMediation_coef_Mcov(X,M.cov,Y,theta)

  return(re)
}
#################################################
