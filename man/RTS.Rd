\name{RTS}
\alias{RTS}
\title{
Robust Two Step Algorithm for Large-Dimensional Elliptical Factor Models
}
\description{
This function is to fit the large-dimensional elliptical factor models via the Robust Two Step (RTS) algorithm.
}
\usage{
RTS(X, r)
}
\arguments{
  \item{X}{
Input matrix, of dimension \eqn{T\times N}. Each row is an observation with \eqn{N} features at time point \eqn{t}.
}
  \item{r}{
A positive integer indicating the factor numbers.
}
}
\details{
See He et al. (2022) for details.
}
\value{
The return value is a list. In this list, it contains the following:
\item{Fhat}{The estimated factor matrix of dimension \eqn{T\times r}.}
\item{Lhat}{The estimated loading matrix of dimension \eqn{N\times r}.}
}
\references{
He, Y., Kong, X., Yu, L., Zhang, X., 2022. Large-dimensional factor analysis without moment constraints.
Journal of Business & Economic Statistics 40, 302–312.
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

fit=RTS(X,3)
fit$Fhat;fit$Lhat
}
