\name{kendall.c}
\alias{kendall.c}
\title{
Calculate modified Kendall's tau distance (native C implementation)
}
\description{
Kendall's tau is equal to the number of adjunct pairwise exchanges required to convert one ranking into another. This modified version allows for partial lists to be compared.
}
\usage{
kendall.c(rank.a, rank.b, k.a, k.b, n, p = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rank.a}{
A single top-\emph{k} list
}
  \item{rank.b}{
A vector of matrix form of top-\emph{k} list(s) to be compared to list \code{a}
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
\value{
Returns modified Kendall's tau distance against \code{a} for each list within \code{b} 
}
\author{
Jie Ding <jding@jimmy.harvard.edu>
}
\note{
An alternative implementation using pure R code is given as \code{\link{kendall}}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{kendall}}, \code{\link{spearman}}
}
\examples{
set.seed(1234)
a <- sample(1:10, 10)
b <- sample(1:10, 10)
kendall.c(a, b, 6, 6, 10)
}

