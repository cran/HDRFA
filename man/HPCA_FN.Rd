\name{HPCA_FN}
\alias{HPCA_FN}
\title{
Estimating Factor Numbers via Rank Minimization Corresponding to Huber PCA
}
\description{
This function is to estimate factor numbers via rank minimization corresponding to Huber Principal Component Analysis (HPCA).
}
\usage{
HPCA_FN(X, rmax, Method = "E", threshold = NULL, L_init = NULL, F_init = NULL, 
        maxiter_HPCA = 100, maxiter_HLM = 100, eps = 0.001)
}
\arguments{
  \item{X}{
Input matrix, of dimension \eqn{T\times N}. Each row is an observation with \eqn{N} features at time point \eqn{t}.
}
  \item{rmax}{
The user-supplied maximum factor numbers. 
}
  \item{Method}{
\code{Method="P"} indicates minimizing the Huber loss of the idiosyncratic error's \eqn{\ell_2} norm while \code{Method="E"} indicates minimizing the elementwise Huber loss. The default is the elementwise Huber loss.
}
  \item{threshold}{
The threshold of rank minimization; default is NULL.
} 
  \item{L_init}{
User-supplied inital value of loadings in the HPCA; default is the PCA estimator.
}
  \item{F_init}{
User-supplied inital value of factors in the HPCA; default is the PCA estimator.
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
\item{rhat}{The estimated factor number.}
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

\donttest{
HPCA_FN(X,8,Method="E")
}
HPCA_FN(X,8,Method="P")
}
