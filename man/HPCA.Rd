\name{HPCA}
\alias{HPCA}
\title{
Huber Principal Component Analysis for Large-Dimensional Factor Models
}
\description{
This function is to fit the factor models via the Huber Principal Component Analysis (HPCA) method. One is based on minimizing the Huber loss of the idiosyncratic error's \eqn{\ell_2}{L2} norm, which turns out to do Principal Component Analysis (PCA) on the weighted sample covariance matrix and thereby named as Huber PCA. The other one is based on minimizing the elementwise Huber loss, which can be solved by an iterative Huber regression algorithm.}

\usage{
HPCA(X, r, Method = "E", tau = NULL, scale_est="MAD", L_init = NULL, 
     F_init = NULL, maxiter_HPCA = 100, maxiter_HLM = 100, eps = 0.001)
}
\arguments{
  \item{X}{
Input matrix, of dimension \eqn{T\times N}. Each row is an observation with \eqn{N} features at time point \eqn{t}.
}
  \item{r}{
A positive integer indicating the factor numbers.
} 
  \item{Method}{
\code{Method="P"} indicates minimizing the Huber loss of the idiosyncratic error's \eqn{\ell_2} norm while \code{Method="E"} indicates minimizing the elementwise Huber loss. The default is the elementwise Huber loss.
}
  \item{tau}{
Optional user-supplied parameter for Huber loss; default is NULL, and \eqn{\tau} is provided by default.
}
  \item{scale_est}{
A parameter for the elementwise Huber loss. \code{scale_est="MAD"} indicates robust variance estimation in each iteration, while \code{scale_est="const"} indicates fixing user-supplied \eqn{\tau}. The default is \code{scale_est="MAD"}.
}
  \item{L_init}{
User-supplied inital value of loadings; default is the PCA estimator.
}
  \item{F_init}{
User-supplied inital value of factors; default is the PCA estimator.
}
  \item{maxiter_HPCA}{
The maximum number of iterations in the HPCA. The default is \eqn{100}.
}
  \item{maxiter_HLM}{
The maximum number of iterations in the iterative Huber regression algorithm. The default is \eqn{100}.
}
  \item{eps}{
The stopping critetion parameter in the HPCA. The default is 1e-3.
}
}
\details{
See He et al. (2023) for details.
}
\value{
The return value is a list. In this list, it contains the following:
\item{Fhat}{The estimated factor matrix of dimension \eqn{T\times r}.}
\item{Lhat}{The estimated loading matrix of dimension \eqn{N\times r}.}
\item{m}{The number of iterations.}
}
\references{
He Y, Li L, Liu D, Zhou W., 2023 Huber Principal Component Analysis for Large-dimensional Factor Models. 
}
\author{
Yong He, Lingxiao Li, Dong Liu, Wenxin Zhou.
}

\examples{
set.seed(1)
T=50;N=50;r=3
L=matrix(rnorm(N*r,0,1),N,r);F=matrix(rnorm(T*r,0,1),T,r)
E=matrix(rnorm(T*N,0,1),T,N)
X=F\%*\%t(L)+E

fit=HPCA(X,r,Method = "E")
fit$Fhat;fit$Lhat

fit=HPCA(X,r,Method = "P")
fit$Fhat;fit$Lhat
}

