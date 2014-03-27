\name{MC.ranks}
\alias{MC.ranks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MC based rank aggregation
}
\description{
Compute aggregate ranks based on the transition matrix from the three Markov Chain algorithms.
}
\usage{
MC.ranks(elements, trans, a, delta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{elements}{
Unique elements of the union of all input lists - second element of the output list from function \code{trans.matrix}
}
  \item{trans}{
One of the three transition matrices build by function \code{trans.matrix} - 4 (5 or 6)
elements of the output list from function \code{trans.matrix}
}
  \item{a}{
Tuning parameter to make sure Markov Chain with the transition matrix is ergodic; parameter value passed from \code{MC}.
}
  \item{delta}{
Convergence criterion for stationary distribution; parameter value passed from \code{MC}.
}
}
\details{
Compute stationary distribution based on a Markov Chain transition matrix built with function \code{trans.matrix}.
}
\value{
A list with 3 components:
\item{comp1}{Number of iterations to reach the stationary distribution}
\item{comp2}{The stationary distribution}
\item{comp3}{The rankings based on the stationary distribution}
}
\references{
Lin, S. (2010)
Space oriented rank-based data integration. Statistical Applications in Genetics and Molecular Biology 9, Article 20.
}
\author{
Shili Lin <shili@stat.osu.edu>
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{MC}, \link{trans.matrix}}
}
