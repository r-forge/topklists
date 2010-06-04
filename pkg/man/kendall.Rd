\name{kendall}
\alias{kendall}
\title{
Calculate modified Kendall's Tau distance
}
\description{
Kendall's Tau is equal to the number of adjunct pairwise exchanges required to convert one ranking into another. This modified version allows for partial lists to be compared.
}
\usage{
kendall(a, b, n = 0, p = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{
a single top k list
}
  \item{b}{
top k list(s) need to be compared to list a
}
  \item{n}{
total number of objects, numbered from 1 to n
}
  \item{p}{
distance added for tied pair (potential problem when p != 0)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
returns modified Kendall's Tau distance against a for each lists within b 
}
\references{
%% ~put references to the literature/web site here ~
%%TODO: ref
}
\author{
Shili Lin <shili@stat.ohio-state.edu>
}
\note{
A faster implementation using native C code is given as \code{\link{kendall.c}}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{kendall.c}}, \code{\link{spearman}}
}
\examples{
set.seed(1234)
rank.pool <- 1:10
a <- sample(rank.pool, 10)
b <- sample(rank.pool, 10)
#for comparing two lists
kendall(a, b)

#for comparing more than two lists make a matrix
c <- sample(rank.pool, 10)
kendall(a, cbind(b, c))
}
\keyword{ arith1 }
