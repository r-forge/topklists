\name{Kendall2Lists}
\alias{Kendall2Lists}
\alias{Kendall2Lists.c}
\title{
Calculate modified Kendall's tau distance
}
\description{
Kendall's tau is equal to the number of adjunct pairwise exchanges required to convert one ranking into another. This modified version allows for partial lists to be compared.
}
\usage{
Kendall2Lists(rank.a, rank.b, k.a, k.b, n, p = 0)
Kendall2Lists.c(rank.a, rank.b, k.a, k.b, n, p = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rank.a}{
A single top-\emph{k} list
}
  \item{rank.b}{
A vector of matrix form of top-\emph{k} list(s) to be compared to the list \code{a}
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
Distance added for tied pair (potential problem when \code{p != 0})
}
}
\details{
There are two implementations available. Pure R code in \code{kendall} and a faster implementation using native C code \code{kendall.c}.
}
\value{
Returns modified Kendall's tau distance against \code{a} for each list within \code{b} 
}
\references{
Lin, S., Ding, J. (2009)
Integration of ranked lists via Cross Entropy Monte Carlo with applications to mRNA and microRNA studies. Biometrics 65, 9-18.
}
\author{
Jie Ding <jding@jimmy.harvard.edu>
}
\seealso{
\code{\link{Spearman}}
}
\examples{
set.seed(1234)
a <- sample(1:10, 10)
b <- sample(1:10, 10)
Kendall2Lists(a, b, 6, 6, 10)
Kendall2Lists.c(a, b, 6, 6, 10)
}
