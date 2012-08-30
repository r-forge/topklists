\name{kendall2.c}
\alias{kendall2.c}
\title{
Modified Kendall's tau Distance (Native C Implementation)
}
\description{
Kendall's tau distance between top K lists with different underlying spaces
}
\usage{
kendall2.c(rank.a, rank.b, k.a, k.b, n, p = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rank.a}{
A single TopKList
}
  \item{rank.b}{
A vector of matrix form of TopKList(s) to be compared to list \code{a}
}
  \item{k.a}{
Value of \code{k} for \code{rank.a}
}
  \item{k.b}{
Value of \code{k} for \code{rank.b}
}
  \item{n}{
Total number of objects, numbered from 1 to n
}
  \item{p}{
Distance added for tied pair (potential problem when p != 0)
}
}
\value{
Returns modified Kendall's Tau distance against a for each lists within b 
}
\author{
Jie Ding <ding@stat.osu.edu>
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{kendall}},\code{\link{kendall.c}}
}
\examples{
set.seed(1234)
a <- sample(1:10, 10)
b <- sample(5:15, 10)
kendall2.c(a, b, 6, 6, 10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
