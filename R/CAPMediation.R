#################################################
# CAP mediation function: finding k components
#' Mediation Analysis with Graph Mediator
#' @description  This function identifies the first \eqn{k} projection directions that satisfies the log-linear model assumption.
#'
#' @param X A \eqn{n \times q} data matrix, the covariate matrix for \eqn{n} subjects with \eqn{q-1} predictors. The first column should be all ones.
#' @param M A list of length \eqn{n}. Each element of the list is a \eqn{T \times p} matrix, representing the data matrix of \eqn{T} observations from \eqn{p} features.
#' @param Y A \eqn{n \times 1} outcome vector, representing the outcome variable for \eqn{n} subjects.
#' @param H A \eqn{p \times p} positive definite matrix used for the theta constraint. Default is \code{NULL}.
#' @param stop.crt A character vector specifying the stopping criteria. Options include \code{"nD"} and \code{"DfD"}.
#' @param nD An integer, specifying the number of directions to be identified. Default is \code{1}.
#' @param DfD.thred A numeric value specifying the threshold for the DfD metric. Default is \code{2}.
#' @param Y.remove A logical value indicating whether to remove identified components from \code{M}. Default is \code{TRUE}.
#' @param max.itr An integer specifying the maximum number of iterations. Default is \code{1000}.
#' @param tol A numeric value specifying the convergence tolerance. Default is \code{1e-4}.
#' @param trace A logical value indicating whether the solution path should be reported. Default is \code{FALSE}.
#' @param score.return A logical value indicating whether the log-variance in the transformed space should be returned. Default is \code{TRUE}.
#' @param theta0.mat A data matrix representing the initial values of \eqn{\theta}. Default is \code{NULL}, in which case the initial values are chosen randomly.
#' @param ninitial An integer specifying the number of different initial values to be tested. If it is greater than 1, multiple initial values will be tested, and the one that yields the minimum objective function will be reported. Default is \code{NULL}.
#' @param seed An integer specifying the random seed for reproducibility. Default is \code{100}.
#' @param verbose A logical value; if \code{TRUE}, the function will print progress updates during execution, indicating the number of components being processed. Useful for monitoring the progress of the algorithm, especially for long computations. Default is \code{TRUE}.

#' @author Yixi Xu, Indiana University School of Medicine, <xuyix@iu.edu>
#'
#' Yi Zhao, Indiana University School of Medicine,<zhaoyi1026@gmail.com>

#' @examples
#' data(env.data.example)
#' X <- get("X", env.data.example)
#' Y <- get("Y", env.data.example)
#' M<-get("M",env.data.example)
#' re<-CAPMediation(X,M,Y,H=NULL)
#'
#' @keywords models
#'

CAPMediation<-function(X,M,Y,H=NULL,stop.crt=c("nD","DfD"),nD=NULL,DfD.thred=2,Y.remove=TRUE,max.itr=1000,tol=1e-4,trace=FALSE,score.return=TRUE,theta0.mat=NULL,ninitial=NULL,seed=100,verbose=TRUE)
{
  # X: (X, W) matrix
  # M: a list of length n, M
  # Y: n by 1 vector
  # H: a p by p positive definite matrix for theta constraint

  if(stop.crt[1]=="nD"&is.null(nD))
  {
    stop.crt<-"DfD"
  }

  n<-length(Y)
  p<-ncol(M[[1]])

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

  #--------------------------------------------
  # First direction
  tm1<-system.time(re1<-CAPMediation_D1_opt(X=X,M=M,Y=Y,H=H,max.itr=max.itr,tol=tol,trace=trace,score.return=score.return,theta0.mat=theta0.mat,ninitial=ninitial,seed=seed))

  Theta.est<-matrix(re1$theta,ncol=1)
  coef.est<-matrix(c(re1$alpha[1],re1$beta,re1$gamma[1],re1$IE,re1$gamma[1]),ncol=1)
  rownames(coef.est)<-c("alpha","beta","gamma","IE","DE")
  if(q>0)
  {
    ancoef.est<-matrix(c(re1$alpha0,re1$alpha[-1],re1$gamma0,re1$gamma[-1]),ncol=1)
    rownames(ancoef.est)<-paste0(rep(c("M","Y"),each=q+1),"-",rep(c("Intercept",X.names[-1]),2))
  }else
  {
    ancoef.est<-matrix(c(re1$alpha0,re1$gamma0),ncol=1)
    rownames(ancoef.est)<-paste0(c("M","Y"),"-","Intercept")
  }
  cp.time<-matrix(as.numeric(tm1[1:3]),ncol=1)
  rownames(cp.time)<-c("user","system","elapsed")

  if(score.return)
  {
    score<-matrix(re1$score,ncol=1)
  }

  if(verbose)
  {
    print(paste0("Component ",ncol(Theta.est)))
  }
  #--------------------------------------------

  if(stop.crt[1]=="nD")
  {
    if(nD>1)
    {
      for(j in 2:nD)
      {
        re.tmp<-NULL
        try(tm.tmp<-system.time(re.tmp<-CAPMediation_Dk(X=X,M=M,Y=Y,H=H,Theta0=Theta.est,Y.remove=Y.remove,max.itr=max.itr,tol=tol,trace=trace,score.return=score.return,
                                                        theta0.mat=theta0.mat,ninitial=ninitial,seed=seed)))

        if(is.null(re.tmp)==FALSE)
        {
          Theta.est<-cbind(Theta.est,re.tmp$theta)
          coef.est<-cbind(coef.est,c(re.tmp$alpha[1],re.tmp$beta,re.tmp$gamma[1],re.tmp$IE,re.tmp$gamma[1]))
          if(q>0)
          {
            ancoef.est<-cbind(ancoef.est,c(re.tmp$alpha0,re.tmp$alpha[-1],re.tmp$gamma0,re.tmp$gamma[-1]))
          }else
          {
            ancoef.est<-cbind(ancoef.est,c(re.tmp$alpha0,re.tmp$gamma0))
          }
          cp.time<-cbind(cp.time,as.numeric(tm.tmp[1:3]))

          if(score.return)
          {
            score<-cbind(score,re.tmp$score)
          }

          if(verbose)
          {
            print(paste0("Component ",ncol(Theta.est)))
          }
        }else
        {
          break
        }
      }
    }

    colnames(Theta.est)=colnames(coef.est)=colnames(ancoef.est)=colnames(cp.time)<-paste0("C",1:ncol(Theta.est))

    cp.time<-cbind(cp.time,apply(cp.time,1,sum))
    colnames(cp.time)[ncol(cp.time)]<-"Total"

    if(score.return)
    {
      colnames(score)<-paste0("C",1:ncol(Theta.est))
    }

    DfD.out<-diag.level(M,Theta.est)
  }
  if(stop.crt[1]=="DfD")
  {
    nD<-1

    DfD.tmp<-1
    while(DfD.tmp<DfD.thred)
    {
      re.tmp<-NULL
      try(tm.tmp<-system.time(re.tmp<-CAPMediation_Dk(X=X,M=M,Y=Y,H=H,Theta0=Theta.est,Y.remove=Y.remove,max.itr=max.itr,tol=tol,trace=trace,score.return=score.return,
                                                      theta0.mat=theta0.mat,ninitial=ninitial,seed=seed)))

      if(is.null(re.tmp)==FALSE)
      {
        nD<-nD+1

        DfD.out<-diag.level(M,cbind(Theta.est,re.tmp$theta))
        DfD.tmp<-DfD.out$avg.level[nD]

        if(DfD.tmp<DfD.thred)
        {
          Theta.est<-cbind(Theta.est,re.tmp$theta)
          coef.est<-cbind(coef.est,c(re.tmp$alpha[1],re.tmp$beta,re.tmp$gamma[1],re.tmp$IE,re.tmp$gamma[1]))
          if(q>0)
          {
            ancoef.est<-cbind(ancoef.est,c(re.tmp$alpha0,re.tmp$alpha[-1],re.tmp$gamma0,re.tmp$gamma[-1]))
          }else
          {
            ancoef.est<-cbind(ancoef.est,c(re.tmp$alpha0,re.tmp$gamma0))
          }
          cp.time<-cbind(cp.time,as.numeric(tm.tmp[1:3]))

          if(score.return)
          {
            score<-cbind(score,re.tmp$score)
          }

          if(verbose)
          {
            print(paste0("Component ",ncol(Theta.est)))
          }
        }
      }else
      {
        break
      }
    }

    colnames(Theta.est)=colnames(coef.est)=colnames(ancoef.est)=colnames(cp.time)<-paste0("C",1:ncol(Theta.est))

    cp.time<-cbind(cp.time,apply(cp.time,1,sum))
    colnames(cp.time)[ncol(cp.time)]<-"Total"

    if(score.return)
    {
      colnames(score)<-paste0("C",1:ncol(Theta.est))
    }

    DfD.out<-diag.level(M,Theta.est)
  }

  theta.orth<-t(Theta.est)%*%Theta.est

  if(score.return)
  {
    re<-list(theta=Theta.est,coef=coef.est,coef.other=ancoef.est,orthogonality=theta.orth,DfD=DfD.out,score=score,time=cp.time)
  }else
  {
    re<-list(theta=Theta.est,coef=coef.est,coef.other=ancoef.est,orthogonality=theta.orth,DfD=DfD.out,time=cp.time)
  }

  return(re)
}
#################################################
