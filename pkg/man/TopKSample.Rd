\name{TopKSample}
\alias{TopKSample}
\alias{TopKSample.c}
\title{
Sampler to generate N top-k lists according to p
}
\description{
Sampler to generate \emph{N} top-\emph{k} lists according to \emph{p}. A function wrapping to a native C implementation is available as well.
}
\usage{
TopKSample(p, N)
TopKSample.c(p, N)
}
\arguments{
  \item{p}{
Matrix of dimension n*(k+1), n is the number of items (to be ranked) and k is the top elements in the joint ranking. Each column is a multinomial probability vector.
}
  \item{N}{
The number of samples
}
}
\details{
A pure R implementation \code{TopKSample} and a native C method \code{TopKSample.c} are available.
}
\value{
N TopKlists
}
\author{
Jie Ding <jding@jimmy.harvard.edu>
}
\note{
By default the C implementation is used due to its better performance.
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
