\name{trans.matrix}
\alias{trans.matrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute transition matrices
}
\description{
Compute transition matrices for all three MC algorithms
}
\usage{
trans.matrix(input, space)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{
%%     ~~Describe \code{input} here~~
}
  \item{space}{
%%     ~~Describe \code{space} here~~
}
}
\details{
Both input and space are lists of the same length=nList
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#get sample data
data(TopKSpaceSampleInput)

trans.matrix(input, space)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
