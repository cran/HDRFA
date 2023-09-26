\name{RTS_FN}
\alias{RTS_FN}
\title{
Estimating Factor Numbers Robustly via Multivariate Kendall’s Tau Eigenvalue Ratios
}
\description{
This function is to estimate factor numbers robustly via multivariate Kendall’s tau eigenvalue ratios.
}
\usage{
RTS_FN(X, rmax)
}
\arguments{
  \item{X}{Input matrix, of dimension \eqn{T\times N}. Each row is an observation with \eqn{N} features at time point \eqn{t}.}
  \item{rmax}{The user-supplied maximum factor numbers.}
}
\details{
See Yu et al. (2019) for details.
}
\value{
\item{rhat}{The estimated factor number.}
}
\references{
Yu, L., He, Y., Zhang, X., 2019. Robust factor number specification for large-dimensional elliptical factor
model. Journal of Multivariate analysis 174, 104543.
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

RTS_FN(X,8)
}

