\name{env.data.example}
\alias{env.data.example}
\docType{data}
\title{
Simulated data
}
\description{ "env.example" is an R environment containing the data generated from the proposed model with \eqn{p=10}.
}
\usage{data("env.data.example")}
\format{
   \describe{
    \item{\code{X}}{a \eqn{n\times q} (X,W) data matrix, the covariate matrix of \eqn{n} subjects with \eqn{q-1} predictors.}
    \item{\code{Y}}{a \eqn{n\times 1} outcome vector, the covariate matrix of \eqn{n} subjects with \eqn{q-1} predictors. }
    \item{\code{M}}{a list of length \eqn{n}. Each list element is a \eqn{T\times p} matrix, the data matrix of \eqn{T} observations from \eqn{p} features.}
    \item{\code{delta}}{a \eqn{n\times p} matrix, the true projection matrix used to generate the data.}
    \item{\code{Sigma}}{a \eqn{p\times p\times n} array, the covariance matrix of the \eqn{n} subjects.}
  }
}
\details{
For subject \eqn{i} observation \eqn{t} (\eqn{i=1,\dots,n}, \eqn{t=1,\dots,T}), \eqn{y_{it}=(y_{it1},\dots,y_{itp})} was generated from a \eqn{p}-dimensional normal distribution with mean zero and covariance \eqn{\Sigma}, where
  \deqn{\Sigma=\Pi\Lambda\Pi,}
  \eqn{\Pi} is an orthonormal matrix and \eqn{\Lambda=\mathrm{diag}(\lambda_{1},\dots,\lambda_{p})} is a diagonal matrix. The eigenvalues \eqn{\lambda_{ij}} (\eqn{j=1,\dots,p}) satisfies the following log-linear model
  \deqn{log(\lambda_{ij})=x_{i}^\top\beta_{j},}
  where \eqn{\beta_{j}} is the \eqn{j}th column of \code{beta}.
}
\examples{
  data(env.data.example)
  X<-get("X",env.data.example)
  Y<-get("Y",env.data.example)
  M<-get("M",env.data.example)
}
\keyword{datasets}

