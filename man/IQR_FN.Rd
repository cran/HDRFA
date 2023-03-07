\name{IQR_FN}
\alias{IQR_FN}
\title{
Estimating Factor Numbers via Rank Minimization Corresponding to IQR
}
\description{
This function is to estimate factor numbers via rank minimization corresponding to Iterative Quantile Regression (IQR).
}
\usage{
IQR_FN(X, rmax, tau, threshold = NULL, max_iter = 100, eps = 10^(-6))
}
\arguments{
  \item{X}{
Input matrix, of dimension \eqn{T\times N}. Each row is an observation with \eqn{N} features at time point \eqn{t}.
}
  \item{rmax}{
The user-supplied maximum factor numbers.
}
  \item{tau}{
The user-supplied quantile level.
}
  \item{threshold}{
The threshold of rank minimization; default is NULL.
}
  \item{max_iter}{
The maximum number of iterations. The default is \eqn{100}.
}
  \item{eps}{
The stopping critetion parameter of the IQR method. The default is 1e-06.
}
}
\details{
See Chen et al. (2021) for more details.
}
\value{
\item{rhat}{The estimated factor number.}
}
\references{
Chen, L., Dolado, J.J., Gonzalo, J., (2021). Quantile factor models. Econometrica 89, 875–910.
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
IQR_FN(X,8,tau)
}

