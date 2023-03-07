\name{IQR}
\alias{IQR}
\title{
Iterative Quantile Regression Methods for Quantile Factor Models
}
\description{
This function is to fit the quantile factor model via the Iterative Quantile Regression (IQR) algorithm.
}
\usage{
IQR(X, r, tau, max_iter = 100, eps = 10^(-6))
}
\arguments{
  \item{X}{
Input matrix, of dimension \eqn{T\times N}. Each row is an observation with \eqn{N} features at time point \eqn{t}.
}
  \item{r}{
A positive integer indicating the factor numbers.
}
  \item{tau}{
The user-supplied quantile level.
}
  \item{max_iter}{
The maximum number of iterations. The default is \eqn{100}.
}
  \item{eps}{
The stopping critetion parameter. The default is 1e-06.
}
}
\details{
See Chen et al. (2021) and He et al. (2023) for details.
}
\value{
The return value is a list. In this list, it contains the following:
\item{Fhat}{The estimated factor matrix of dimension \eqn{T\times r}.}
\item{Lhat}{The estimated loading matrix of dimension \eqn{N\times r}.}
\item{t}{The number of iterations.}
}
\references{
Chen, L., Dolado, J.J., Gonzalo, J., (2021). Quantile factor models. Econometrica 89, 875–910.

He Y, Kong X, Yu L, Zhao P., (2023). Quantile factor analysis for large-dimensional time series with statistical guarantee <arXiv:2006.08214>.
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

tau=0.5
fit=IQR(X,r,tau)
fit$Fhat;fit$Lhat
}
