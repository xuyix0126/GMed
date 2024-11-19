#################################################
# Bootstrap inference
#' Inference of model coefficients
#' @description This function performs inference on the model coefficient \eqn{\theta}.
#'
#' @param X A \eqn{n \times q} data matrix, the covariate matrix for \eqn{n} subjects with \eqn{q-1} predictors. The first column should be all ones.
#' @param M A list of length \eqn{n}. Each element of the list is a \eqn{T \times p} matrix, representing the data matrix of \eqn{T} observations from \eqn{p} features.
#' @param Y A \eqn{n \times 1} outcome vector, representing the outcome variable for \eqn{n} subjects.
#' @param theta a \eqn{p}-dimensional vector, the projecting direction \eqn{\theta}. Default is \code{NULL}. If \code{theta = NULL}, an error warning will be returned.
#' @param H A \eqn{p \times p} positive definite matrix used for the theta constraint. Default is \code{NULL}.
#' @param boot  a logic variable, whether bootstrap inference is performed.
#' @param sims a numeric value, the number of bootstrap iterations will be performed.
#' @param boot.ci.type a character of the way of calculating bootstrap confidence interval. If \code{boot.ci.type = "bca"}, the bias corrected confidence interval is returned; if \code{boot.ci.type = "perc"}, the percentile confidence interval is returned.
#' @param conf.level a numeric value, the designated significance level. Default is \eqn{0.95}.
#' @param seed.boot An integer specifying the random seed for reproducibility. Default is \code{100}.
#' @param verbose a logic variable, whether the bootstrap procedure is printed. Default is \code{TRUE}.

#' @author Yixi Xu, Indiana University School of Medicine, <xuyix@iu.edu>
#'
#' Yi Zhao, Indiana University School of Medicine,<zhaoyi1026@gmail.com>
#'
#'
#' @examples
#' data(env.data.example)
#' X <- get("X", env.data.example)
#' Y <- get("Y", env.data.example)
#' M<-get("M",env.data.example)
#'
#' gamma.mat0<-matrix(runif(p),nrow=p,ncol=p)
#' gamma.mat<-qr.Q(qr(gamma.mat0))
#' for(j in 1:p)
#' {
#' if(gamma.mat[which.max(abs(gamma.mat[,j])),j]<0)
#' {
#'  gamma.mat[,j]<-(-gamma.mat[,j])
#'  }
#'  }
#'  Gamma<-gamma.mat
#' # get bootstrap result
#' # re.boot<-CAPMediation_boot(X,M,Y,theta=Gamma[,2],H=NULL)
#'
#'
#' @keywords models
#'
#'

CAPMediation_boot<-function(X,M,Y,theta=NULL,H=NULL,boot=TRUE,sims=1000,boot.ci.type=c("se","perc"),conf.level=0.95,seed.boot=100,verbose=TRUE)
{
  # X: (X, W) matrix
  # M: a list of length n, M
  # Y: n by 1 vector
  # theta: components
  # H: a p by p positive definite matrix for theta constraint

  n<-length(Y)
  p<-ncol(M[[1]])

  nX<-ncol(X)
  q<-ncol(X)-1
  if(is.null(colnames(X)))
  {
    if(q>0)
    {
      X.names<-c("X",paste0("W",1:q))
    }else
    {
      X.names<-c("X")
    }
    colnames(X)<-X.names
  }else
  {
    X.names<-colnames(X)
  }

  M.cov<-array(NA,c(p,p,n))
  for(i in 1:n)
  {
    M.cov[,,i]<-cov(M[[i]])
  }

  if(is.null(H))
  {
    H<-apply(M.cov,c(1,2),mean,na.rm=TRUE)
  }

  if(boot)
  {
    if(is.null(theta))
    {
      stop("Error! Need theta value.")
    }else
    {
      coef.boot<-matrix(NA,5,sims)
      rownames(coef.boot)<-c("alpha","beta","gamma","IE","DE")
      if(q>0)
      {
        coef.other.boot<-matrix(NA,(q+1)*2,sims)
        rownames(coef.other.boot)<-paste0(rep(c("M","Y"),each=q+1),"-",rep(c("Intercept",X.names[-1]),2))
      }else
      {
        coef.other.boot<-matrix(NA,2,sims)
        rownames(coef.other.boot)<-paste0(c("M","Y"),"-","Intercept")
      }

      for(b in 1:sims)
      {
        set.seed(seed.boot+b)
        idx.tmp<-sample(1:n,n,replace=TRUE)

        Ytmp<-Y[idx.tmp]
        Xtmp<-matrix(X[idx.tmp,],ncol=nX)
        colnames(Xtmp)<-X.names
        Mtmp<-M[idx.tmp]

        re.tmp<-NULL
        try(re.tmp<-CAPMediation_coef(X=Xtmp,M=Mtmp,Y=Ytmp,theta=theta))
        if(is.null(re.tmp)==FALSE)
        {
          coef.boot[,b]<-c(re.tmp$alpha[1],re.tmp$beta,re.tmp$gamma[1],re.tmp$IE,re.tmp$gamma[1])
          coef.other.boot[,b]<-c(re.tmp$alpha0,re.tmp$alpha[-1],re.tmp$gamma0,re.tmp$gamma[-1])
        }

        if(verbose)
        {
          print(paste0("Bootstrap sample ",b))
        }
      }

      coef.est<-matrix(NA,nrow(coef.boot),6)
      rownames(coef.est)<-rownames(coef.boot)
      colnames(coef.est)<-c("Estimate","SE","statistics","pvalue","LB","UB")
      coef.est[,1]<-apply(coef.boot,1,mean,na.rm=TRUE)
      coef.est[,2]<-apply(coef.boot,1,sd,na.rm=TRUE)
      coef.est[,3]<-coef.est[,1]/coef.est[,2]
      coef.est[,4]<-(1-pnorm(abs(coef.est[,3])))*2

      coef.other.est<-matrix(NA,nrow(coef.other.boot),6)
      rownames(coef.other.est)<-rownames(coef.other.boot)
      colnames(coef.other.est)<-c("Estimate","SE","statistics","pvalue","LB","UB")
      coef.other.est[,1]<-apply(coef.other.boot,1,mean,na.rm=TRUE)
      coef.other.est[,2]<-apply(coef.other.boot,1,sd,na.rm=TRUE)
      coef.other.est[,3]<-coef.other.est[,1]/coef.other.est[,2]
      coef.other.est[,4]<-(1-pnorm(abs(coef.other.est[,3])))*2

      if(boot.ci.type[1]=="se")
      {
        coef.est[,5]<-coef.est[,1]-qnorm(1-(1-conf.level)/2)*coef.est[,2]
        coef.est[,6]<-coef.est[,1]+qnorm(1-(1-conf.level)/2)*coef.est[,2]

        coef.other.est[,5]<-coef.other.est[,1]-qnorm(1-(1-conf.level)/2)*coef.other.est[,2]
        coef.other.est[,6]<-coef.other.est[,1]+qnorm(1-(1-conf.level)/2)*coef.other.est[,2]
      }
      if(boot.ci.type[1]=="perc")
      {
        coef.est[,5]<-apply(coef.boot,1,quantile,probs=(1-conf.level)/2)
        coef.est[,6]<-apply(coef.boot,1,quantile,probs=1-(1-conf.level)/2)

        coef.other.est[,5]<-apply(coef.other.boot,1,quantile,probs=(1-conf.level)/2)
        coef.other.est[,6]<-apply(coef.other.boot,1,quantile,probs=1-(1-conf.level)/2)
      }
    }

    re<-list(coef=coef.est,coef.other=coef.other.est,coef.boot=coef.boot,coef.other.boot=coef.other.boot)

    return(re)
  }else
  {
    stop("Error!")
  }
}
#################################################

#################################################
# CAP mediation refit with k components
CAPMediation_refit_Mcov<-function(X,M.cov,Y,Theta)
{
  # X: (X, W) matrix
  # M.cov: p by p by n array, covariance matrix of M
  # Y: n by 1 vector
  # Theta: p by k projection matrix

  k<-ncol(Theta)
  if(k==1)
  {
    re<-CAPMediation_coef_Mcov(X,M.cov,Y,theta=Theta)
  }else
  {
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

    score<-matrix(NA,n,k)
    colnames(score)<-paste0("C",1:k)
    for(jj in 1:k)
    {
      score[,jj]<-apply(M.cov,c(3),function(x){return((t(Theta[,jj])%*%x%*%Theta[,jj])[1,1])})
    }

    # beta and gamma estimate
    Z<-cbind(rep(1,n),X,log(score))
    colnames(Z)<-c("Intercept",colnames(X),paste0("M",1:k))
    mu.est<-c(ginv(t(Z)%*%Z)%*%t(Z)%*%Y)

    gamma0.est<-mu.est[1]
    gamma.est<-mu.est[2:(length(mu.est)-k)]
    beta.est<-mu.est[(length(mu.est)-k+1):length(mu.est)]

    # sigma2 estimate
    sigma2.est<-mean((Y-Z%*%mu.est)^2,na.rm=TRUE)

    # estimate alpha0, alpha0.rnd, alpha using mixed effects model
    alpha0.est<-rep(NA,k)
    alpha.est<-matrix(NA,ncol(X),k)
    rownames(alpha.est)<-colnames(X)
    colnames(alpha.est)<-paste0("M",1:k)
    alpha0.rnd.est<-matrix(NA,n,k)
    colnames(alpha0.rnd.est)<-paste0("M",1:k)
    tau2.est<-rep(NA,k)
    for(jj in 1:k)
    {
      dtmp<-data.frame(ID=1:n,score=log(score[,jj]),X)
      eval(parse(text=paste0("fit.tmp<-lme(score~",paste(colnames(X),collapse="+"),",random=~1|ID,data=dtmp,control=lmeControl(opt='optim'))")))

      alpha0.est[jj]<-fit.tmp$coefficients$fixed[1]
      alpha.est[,jj]<-fit.tmp$coefficients$fixed[-1]
      alpha0.rnd.est[,jj]<-c(fit.tmp$coefficients$random$ID+fit.tmp$coefficients$fixed[1])
      tau2.est[jj]<-mean((alpha0.rnd.est[,jj]-alpha0.est[jj])^2,na.rm=TRUE)
    }

    IE.est<-alpha.est[1,]*beta.est

    names(alpha0.est)=names(beta.est)=names(tau2.est)<-paste0("M",1:k)
    names(gamma.est)<-colnames(X)
    re<-list(theta=Theta,alpha=alpha.est,beta=beta.est,gamma=gamma.est,IE=IE.est,alpha0=alpha0.est,alpha0.rnd=alpha0.rnd.est,gamma0=gamma0.est,tau2=tau2.est,sigma2=sigma2.est)

    return(re)
  }
}
CAPMediation_refit<-function(X,M,Y,Theta)
{
  # X: (X, W) matrix
  # M: a list of length n, M
  # Y: n by 1 vector
  # Theta: p by k projection matrix

  M.cov<-array(NA,c(p,p,n))
  for(i in 1:n)
  {
    M.cov[,,i]<-cov(M[[i]])
  }

  re<-CAPMediation_refit_Mcov(X,M.cov,Y,Theta)

  return(re)
}
#################################################

#################################################
# CAP mediation refit bootstrap inference
CAPMediation_refit_boot<-function(X,M,Y,Theta=NULL,H=NULL,boot=TRUE,sims=1000,boot.ci.type=c("se","perc"),conf.level=0.95,seed.boot=100,verbose=TRUE)
{
  # X: (X, W) matrix
  # M: a list of length n, M
  # Y: n by 1 vector
  # Theta: p by k projection matrix
  # H: a p by p positive definite matrix for theta constraint

  if(is.null(Theta)|ncol(Theta)==1)
  {
    re<-CAPMediation_boot(X,M,Y,theta=Theta,H=H,boot=boot,sims=sims,boot.ci.type=boot.ci.type,conf.level=conf.level,seed.boot=seed.boot,verbose=verbose)
  }else
  {
    k<-ncol(Theta)

    n<-length(Y)
    p<-ncol(M[[1]])

    nX<-ncol(X)
    q<-ncol(X)-1
    if(is.null(colnames(X)))
    {
      if(q>0)
      {
        X.names<-c("X",paste0("W",1:q))
      }else
      {
        X.names<-c("X")
      }
      colnames(X)<-X.names
    }else
    {
      X.names<-colnames(X)
    }

    M.cov<-array(NA,c(p,p,n))
    for(i in 1:n)
    {
      M.cov[,,i]<-cov(M[[i]])
    }

    if(is.null(H))
    {
      H<-apply(M.cov,c(1,2),mean,na.rm=TRUE)
    }

    if(boot)
    {
      coef.boot<-matrix(NA,3*k+2,sims)
      rownames(coef.boot)<-c(paste0("alpha_M",1:k),paste0("beta_M",1:k),"gamma",paste0("IE_M",1:k),"DE")
      if(q>0)
      {
        coef.other.boot<-matrix(NA,(q+1)*(k+1),sims)
        rownames(coef.other.boot)<-paste0(rep(c(paste0("M",1:k),"Y"),each=q+1),"-",rep(c("Intercept",X.names[-1]),k+1))
      }else
      {
        coef.other.boot<-matrix(NA,k+1,sims)
        rownames(coef.other.boot)<-paste0(c(paste0("M",1:k),"Y"),"-","Intercept")
      }

      for(b in 1:sims)
      {
        set.seed(seed.boot+b)
        idx.tmp<-sample(1:n,n,replace=TRUE)

        Ytmp<-Y[idx.tmp]
        Xtmp<-matrix(X[idx.tmp,],ncol=nX)
        colnames(Xtmp)<-X.names
        Mtmp<-M[idx.tmp]

        re.tmp<-NULL
        try(re.tmp<-CAPMediation_refit(X=Xtmp,M=Mtmp,Y=Ytmp,Theta=Theta))
        if(is.null(re.tmp)==FALSE)
        {
          coef.boot[,b]<-c(re.tmp$alpha[1,],re.tmp$beta,re.tmp$gamma[1],re.tmp$IE,re.tmp$gamma[1])
          coef.other.boot[,b]<-c(c(rbind(re.tmp$alpha0,re.tmp$alpha[-1,])),re.tmp$gamma0,re.tmp$gamma[-1])
        }

        if(verbose)
        {
          print(paste0("Bootstrap sample ",b))
        }
      }

      coef.est<-matrix(NA,nrow(coef.boot),6)
      rownames(coef.est)<-rownames(coef.boot)
      colnames(coef.est)<-c("Estimate","SE","statistics","pvalue","LB","UB")
      coef.est[,1]<-apply(coef.boot,1,mean,na.rm=TRUE)
      coef.est[,2]<-apply(coef.boot,1,sd,na.rm=TRUE)
      coef.est[,3]<-coef.est[,1]/coef.est[,2]
      coef.est[,4]<-(1-pnorm(abs(coef.est[,3])))*2

      coef.other.est<-matrix(NA,nrow(coef.other.boot),6)
      rownames(coef.other.est)<-rownames(coef.other.boot)
      colnames(coef.other.est)<-c("Estimate","SE","statistics","pvalue","LB","UB")
      coef.other.est[,1]<-apply(coef.other.boot,1,mean,na.rm=TRUE)
      coef.other.est[,2]<-apply(coef.other.boot,1,sd,na.rm=TRUE)
      coef.other.est[,3]<-coef.other.est[,1]/coef.other.est[,2]
      coef.other.est[,4]<-(1-pnorm(abs(coef.other.est[,3])))*2

      if(boot.ci.type[1]=="se")
      {
        coef.est[,5]<-coef.est[,1]-qnorm(1-(1-conf.level)/2)*coef.est[,2]
        coef.est[,6]<-coef.est[,1]+qnorm(1-(1-conf.level)/2)*coef.est[,2]

        coef.other.est[,5]<-coef.other.est[,1]-qnorm(1-(1-conf.level)/2)*coef.other.est[,2]
        coef.other.est[,6]<-coef.other.est[,1]+qnorm(1-(1-conf.level)/2)*coef.other.est[,2]
      }
      if(boot.ci.type[1]=="perc")
      {
        coef.est[,5]<-apply(coef.boot,1,quantile,probs=(1-conf.level)/2)
        coef.est[,6]<-apply(coef.boot,1,quantile,probs=1-(1-conf.level)/2)

        coef.other.est[,5]<-apply(coef.other.boot,1,quantile,probs=(1-conf.level)/2)
        coef.other.est[,6]<-apply(coef.other.boot,1,quantile,probs=1-(1-conf.level)/2)
      }

      re<-list(coef=coef.est,coef.other=coef.other.est,coef.boot=coef.boot,coef.other.boot=coef.other.boot)

      return(re)
    }else
    {
      stop("Error!")
    }
  }
}
#################################################

