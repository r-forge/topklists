\name{kendall.c}
\alias{kendall.c}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate modified Kendall's Tau (using native C code)
}
\description{
Kendall's Tau is equal to the number of adjunct pairwise exchanges required to convert one ranking into another. This modified version allows for partial lists to be compared.
}
\usage{
kendall.c(a, b, n = 0, p = 0)
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
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
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
\keyword{ ~kwd1 }

