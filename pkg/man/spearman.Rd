\name{spearman}
\alias{spearman}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate modified Spearman's footrule distance between two top K lists
}
\description{
Spearman's footrule is a measure for distance between ranked lists. It is given as the sum of absolute differences between ranks of two lists. This modificated version allows for comparing partial lists.
}
\usage{
spearman(rank.a, rank.b, k.a, k.b, n)
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
}
\value{
Returns modified spearman distance against \code{a} for each lists within \code{b} 
}
\author{
Jie Ding <jding@jimmy.harvard.edu>
}
\note{
Implemented using native C code. For the R version see \code{\link{kendall}}.
}
\seealso{
\code{\link{kendall}}, \code{\link{spearman}}
}
\examples{
set.seed(1234)
a <- sample(1:10, 10)
b <- sample(1:10, 10)
spearman(a, b, 6, 6, 10)

}
\keyword{ arith }

