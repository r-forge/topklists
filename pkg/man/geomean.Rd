\name{geomean}
\alias{geomean}
\title{
Calculate the geometric mean
}
\description{
Calculate the geometric mean
}
\usage{
geomean(x, na.rm = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a vector of values
}
  \item{na.rm}{
Should missing values automatically be removed fro calculation
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
the geometric mean
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Shili Lin <shili@stat.ohio-state.edu>
}
\note{
%%  ~~further notes~~
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
set.seed(122)
vals <- sample(1:100, 10)
geomean(vals)
}
\keyword{ arith }

