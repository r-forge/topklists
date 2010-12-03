\name{prepareIdata}
\alias{prepareIdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prepare data from rankings
}
\description{
function to prepare Idata from several rankings of several assessors
}
\usage{
prepareIdata(x, delta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Data matrix, where columns represents rank order of genes from two different assessors, rows represents genes (clones)
}
  \item{delta}{
The maximal distance of the objects in ranking between assessors}
}
\details{
}
\value{
The result is a object of type "Idata", which is a list containing Idata and the information
about the distance (delta).
\item{Idata}{ranking distance information}
\item{delta}{the given delta}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Eva Budinska <budinska@iba.muni.cz>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
set.seed(4568)
A <- sample(1:20, 20)
B <- sample(1:20, 20)
C <- sample(1:20, 20)
mm <- data.frame(A, B, C, row.names=LETTERS[1:20])
prepareIdata(mm, delta=10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ list }

