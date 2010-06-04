\name{comp.rank}
\alias{comp.rank}
\title{
Calculate composite rank
}
\description{
Calculate composite rank from an input list of top-k lists.
}
\usage{
comp.rank(input.list)
}
\arguments{
  \item{input.list}{
a list of topk lists, may have different length
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns the composite ranks built from input lists
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Shili Lin <shili@stat.ohio-state.edu>
}
\note{
}
\seealso{
}
\examples{
set.seed(1234)

rank.pool <- 1:10
a <- sample(rank.pool, 10)
b <- sample(rank.pool, 10)
c <- sample(rank.pool, 10)

comp.rank(list(a, b, c))
}
\keyword{ robust }
