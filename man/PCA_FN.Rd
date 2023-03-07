\name{PCA_FN}
\alias{PCA_FN}
\title{
Estimating Factor Numbers via Eigenvalue Ratios Corresponding to PCA
}
\description{
This function is to estimate factor numbers via eigenvalue ratios corresponding to Principal Component Analysis (PCA).
}
\usage{
PCA_FN(X, rmax)
}
\arguments{
  \item{X}{Input matrix, of dimension \eqn{T\times N}. Each row is an observation with \eqn{N} features at time point \eqn{t}.}
  \item{rmax}{The user-supplied maximum factor numbers.}
}
\details{
See Ahn and Horenstein (2013) for details.
}
\value{
\item{rhat}{The estimated factor numbers.}

}
\references{
Ahn, S.C., Horenstein, A.R., (2013). Eigenvalue ratio test for the number of factors. Econometrica 81,
1203–1227.
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

PCA_FN(X,8)
}
