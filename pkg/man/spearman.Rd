\name{spearman}
\alias{spearman}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate modified Spearman's footrule distance between two top K lists
}
\description{
Spearman's footrule is a measure fir distance between ranked lists. It given as the sum of absolute differences between ranks of two lists. This modificated version allows for comparing partial lists.
}
\usage{
spearman(a, b, n = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{
a single top k list
}
  \item{b}{
a matrix with the top k list(s) need to be compared to list a, each column is one list
}
  \item{n}{
total number of objects, numbered from 1 to n
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
returns modified spearman distance against a for each lists within b 
}
\references{
%% ~put references to the literature/web site here ~
%%TODO: ref
}
\author{
Shili Lin <shili@stat.ohio-state.edu>
}
\note{
Implemented using native C code. For the R version see \code{\link{kendall}}.
}
\seealso{
\code{\link{kendall}}, \code{\link{spearman}}
}
\examples{
set.seed(1234)
rank.pool <- 1:10
a <- sample(rank.pool, 10)
b <- sample(rank.pool, 10)
#for comparing two lists
kendall.c(a, b)

#for comparing more than two lists make a matrix
c <- sample(rank.pool, 10)
kendall.c(a, cbind(b, c))

}
\keyword{ arith }

