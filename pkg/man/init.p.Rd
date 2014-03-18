\name{init.p}
\alias{init.p}
\title{
Initialization method for probabilities
}
\description{
Initialization method for probabilities
}
\usage{
init.p(topK, n, k, init.m = "p", init.w = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{topK}{
A list of input lists, with items coded from 1 to \code{n}
}
  \item{n}{
Total number of items
}
  \item{k}{
Length of target list
}
  \item{init.m}{
Initialization method, "p" point mass, "s" smooth, "cp" point mass using composite ranks, "cs" smooth using composite ranks
}
  \item{init.w}{
initialization weight
}
}
\value{
A probability matrix
}
\references{
Lin, S., Ding, J. (2009)
Integration of Ranked Lists  via
Cross Entropy Monte Carlo with Applications to mRNA and microRNA
Studies. Biometrics 65, 9-18.
}
\author{
Jie Ding <jding@jimmy.harvard.edu>
}
\examples{
set.seed(1234)

rank.pool <- 1:10
a <- sample(rank.pool, 10)
b <- sample(rank.pool, 10)
c <- sample(rank.pool, 10)
rlist <- list(a, b, c)


init.p(rlist, length(unique(unlist(rlist))), 5, "cp")
}
