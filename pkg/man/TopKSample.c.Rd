\name{TopKSample.c}
\alias{TopKSample.c}
\title{
Sampler to generate N top K lists according to p (implemented in native C)
}
\description{
Sampler to generate N top K lists according to p. This function is implemented using C. A pure R version is available as well.
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
N TopKlists}
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
