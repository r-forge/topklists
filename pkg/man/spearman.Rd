\name{Spearman}
\alias{Spearman}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Modified Spearman's footrule distance
}
\description{
Spearman's footrule is a measure for distance between ranked lists. It is given as the sum of absolute differences between ranks of two lists. Here a modified version is implemented that allows for comparing partial lists.
}
\usage{
Spearman(rank.a, rank.b, k.a, k.b, n)
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
}
\value{
Returns modified Spearman distance against \code{a} for each lists within \code{b} 
}
\author{
Jie Ding <jding@jimmy.harvard.edu>
}
\seealso{
\code{\link{Kendall2Lists}}
}
\examples{
set.seed(1234)
a <- sample(1:10, 10)
b <- sample(1:10, 10)
Spearman(a, b, 6, 6, 10)

}

