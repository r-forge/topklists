\name{kendall}
\alias{kendall}
\title{
Calculate modified Kendall's Tau distance
}
\description{
Kendall's Tau is equal to the number of adjunct pairwise exchanges required to convert one ranking into another. This modified version allows for partial lists to be compared.
}
\usage{
kendall(rank.a, rank.b, k.a, k.b, n, p = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rank.a}{
A single top-K list
}
  \item{rank.b}{
A vector of matrix form of top-K list(s) to be compared to list \code{a}
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
Returns modified Kendall's Tau distance against \code{a} for each lists within \code{b} 
}
\references{
Lin, S., Ding, J. (2009)
Integration of ranked lists via Cross Entropy Monte Carlo with applications to mRNA and microRNA studies. Biometrics 65, 9-18.
}
\author{
Jie Ding <jding@jimmy.harvard.edu>
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
a <- sample(1:10, 10)
b <- sample(1:10, 10)
kendall(a, b, 6, 6, 10)
}
