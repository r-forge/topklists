\name{TopKSample.c}
\alias{TopKSample.c}
\title{
Sampler to generate N top-k lists according to p (implemented in native C)
}
\description{
Sampler to generate \emph{N} top-\emph{k} lists according to \code{p}. This function is implemented using C. A pure R version is available as well.
}
\usage{
TopKSample.c(p, N)
}
\arguments{
  \item{p}{
Matrix of dimension n*(k+1), n is the number of items
}
  \item{N}{
The number of samples
}
}
\value{
N top-K lists}
\author{
Jie Ding <jding@jimmy.harvard.edu>
}
\note{
By default the C implementation is used due to its better performance.
}
\seealso{
\code{\link{TopKSample}}
}
\examples{
set.seed(1234)
rank.pool <- 1:10
a <- sample(rank.pool, 10)
b <- sample(rank.pool, 10)
c <- sample(rank.pool, 10)
M <- cbind(a, b, c)

TopKSample.c(M, 4)

}

