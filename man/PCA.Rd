\name{PCA}
\alias{PCA}
\title{
Principal Component Analysis for Large-Dimensional Factor Models
}
\description{
This function is to fit the factor models via Principal Component Analysis (PCA) methods.
}
\usage{
PCA(X, r, constraint = "L")
}
\arguments{
  \item{X}{
Input matrix, of dimension \eqn{T\times N}. Each row is an observation with \eqn{N} features at time point \eqn{t}.
}
  \item{r}{
A positive integer indicating the factor numbers.
}
  \item{constraint}{
The type of identification condition. If \code{constraint="L"}, the columns of the estimated loading matrix are orthogonal and \code{constraint="F"} indicates the columns of the estimated factor matrix are orthogonal.
}
}
\details{
See Bai (2003) for details.
}
\value{
The return value is a list. In this list, it contains the following:
\item{Fhat}{The estimated factor matrix of dimension \eqn{T\times r}.}
\item{Lhat}{The estimated loading matrix of dimension \eqn{N\times r}.}
}
\references{
Bai, J., 2003. Inferential theory for factor models of large dimensions. Econometrica 71, 135–171.
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

fit=PCA(X,3,"L")
t(fit$Lhat)\%*\%fit$Lhat/N

fit=PCA(X,3,"F")
t(fit$Fhat)\%*\%fit$Fhat/T
}

