#################################################
# eigenvectors and eigenvalues of A with respect to H
# H positive definite and symmetric
eigen.solve<-function(A,H)
{
  p<-ncol(H)

  H.svd<-svd(H)
  H.d.sqrt<-diag(sqrt(H.svd$d))
  H.d.sqrt.inv<-diag(1/sqrt(H.svd$d))
  H.sqrt.inv<-H.svd$u%*%H.d.sqrt.inv%*%t(H.svd$u)

  #---------------------------------------------------
  # svd decomposition method
  eigen.tmp<-eigen(H.d.sqrt.inv%*%t(H.svd$u)%*%A%*%H.svd$u%*%H.d.sqrt.inv)
  eigen.tmp.vec<-Re(eigen.tmp$vectors)

  # obj<-rep(NA,ncol(eigen.tmp$vectors))
  # for(j in 1:ncol(eigen.tmp$vectors))
  # {
  #   otmp<-H.svd$u%*%H.d.sqrt.inv%*%eigen.tmp.vec[,j]
  #   obj[j]<-t(otmp)%*%A%*%otmp
  # }
  # re<-H.svd$u%*%H.d.sqrt.inv%*%eigen.tmp.vec[,which.min(obj)]
  re<-H.svd$u%*%H.d.sqrt.inv%*%eigen.tmp.vec[,p]
  #---------------------------------------------------

  #---------------------------------------------------
  # eigenvector of A with respect to H
  # eigen.tmp<-eigen(H.sqrt.inv%*%A%*%H.sqrt.inv)
  # eigen.tmp.vec<-Re(eigen.tmp$vectors)
  #
  # obj<-rep(NA,ncol(eigen.tmp$vectors))
  # for(j in 1:ncol(eigen.tmp$vectors))
  # {
  #   otmp<-H.sqrt.inv%*%eigen.tmp.vec[,j]
  #   obj[j]<-t(otmp)%*%A%*%otmp
  # }
  # re<-H.sqrt.inv%*%eigen.tmp.vec[,p]
  #---------------------------------------------------

  return(c(re))
}
#################################################
