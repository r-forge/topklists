\name{compute.stream}
\alias{compute.stream}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The actual function estimating the point of degeneration on the signal-noise vector of 0's and 1's created by prepareIdata.
}
\description{
The algorithm consists of an ordered sequence of "test stages" s1, s2....
I stage sk an integer Jsk is estimated, which is a potential lower bound to j0 (when k is odd), or potential upper bound to j0 (when k is even).
}
\usage{
compute.stream(Idata, const=0.251, v, r=1.2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Idata}{
Input data is a signal-noise vector of 0s and 1s (see \code{prepareIdata})
}
  \item{const}{
const is 
}
  \item{v}{
Represents a tuning parameter representing the degree of randomness in the assignments. In each step the noise is estimated from Idata as probability of 1 within the interval of size v, moving from Jsk\-v\*r if k is odd or Jsk\+v\*r if k is even, until the convergence or a break (see r).
}
  \item{r}{
a constant determining the starting point from which the probability of 1 is estimated in a winodow of size v (see \code{v}).
}
}
\value{
A named list containing:
 \item{j0_est }{The estimated position for Idata turning into noise}
 \item{reason.break}{The reason why the computation has ended - convergence or other break}
 \item{js}{The sequence of estimated j0 in each round, shows also the convergence}
 \item{v}{Value of parameter v used for computation}
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

