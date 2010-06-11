\name{j0.multi}
\alias{j0.multi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Method for counting K(j0_est)
}
\description{
Method to count K(j0_est) for each pair of several ranked lists. 
}
\usage{
j0.multi(Data, m.dist, v)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{
input matrix, each coloum presents one list
}
  \item{m.dist}{
the maximal distance of the objects in ranking
}
  \item{v}{
parameter for estimating j0
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a list containing
\item{maxK}{maximal estimated J0 through all combinations of 2 lists}
\item{K}{matrix with all pair-wise comparison}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
set.seed(4657)

A <- sample(1:10, 10)
B <- sample(1:10, 10)
C <- sample(1:10, 10)

j0.multi(cbind(A, B, C), 10, 10)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ list }

