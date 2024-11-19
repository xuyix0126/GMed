#################################################
# level of diagonalization
diag.level<-function(M,Theta)
{
  # M: a list of length n, M
  # Theta: p by k matrix, identified components

  n<-length(M)

  if(ncol(Theta)==1)
  {
    re<-list(avg.level=1,sub.level=rep(1,n))
  }else
  {
    p<-ncol(M[[1]])
    nT<-sapply(M,nrow)
    ps<-ncol(Theta)

    dl.sub<-matrix(NA,n,ps)
    colnames(dl.sub)<-paste0("C",1:ps)
    dl.sub[,1]<-1
    for(i in 1:n)
    {
      cov.tmp<-cov(M[[i]])

      for(j in 2:ps)
      {
        theta.tmp<-Theta[,1:j]
        mat.tmp<-t(theta.tmp)%*%cov.tmp%*%theta.tmp
        dl.sub[i,j]<-det(diag(diag(mat.tmp)))/det(mat.tmp)
      }
    }

    pmean<-apply(dl.sub,2,function(y){return(prod(apply(cbind(y,nT),1,function(x){return(x[1]^(x[2]/sum(nT)))})))})

    re<-list(avg.level=pmean,sub.level=dl.sub)
  }

  return(re)
}
#################################################
