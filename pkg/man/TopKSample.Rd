\name{TopKSample}
\alias{TopKSample}
\title{
Sampler to generate N top K lists according to p (implemented in native C)
}
\description{
Sampler to generate N top K lists according to p. A function wrapping to a native C implementation is available as well.
}
\usage{
TopKSample(p, N)
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
N TopKlists
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Jie Ding <ding@stat.osu.edu>
}
\note{
By default the C implementation is used due to the better performance.
}
\seealso{
\code{\link{TopKSample.c}}
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
\keyword{ distribution }
