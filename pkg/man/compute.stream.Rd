\name{compute.stream}
\alias{compute.stream}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Method to estimate j0
}
\description{
The algorithm consists of an ordered sequence of "test stages" s1, s2....
Stage sk is associated with an integer Jsk, which when k is odd, is a potential lower bound to j0
}
\usage{
compute.stream(Idata, const=0.251, v, r=1.2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Idata}{
Input data is a vector of 0's and 1's
}
  \item{const}{
%%     ~~Describe \code{const} here~~
}
  \item{v}{
%%     ~~Describe \code{v} here~~
}
  \item{r}{
%%     ~~Describe \code{r} here~~
}
}
\value{
A named list containing:
 \item{j0_est }{The estimated position for Idata turning into noise}
 \item{reason.break}{The reason why the computation has ended}
 \item{js}{}
 \item{v.vector}{}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
 Eva Budinska <budinska@iba.muni.cz>
}
\examples{
set.seed(465)
myhead <- rbinom(20, 1, 0.8)
mytail <- rbinom(20, 1, 0.5)
mydata <- c(myhead, mytail)
compute.stream(mydata, v=10)	
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ list}

