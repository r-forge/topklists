\name{MC.ranks}
\alias{MC.ranks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute rankings based on a MC algorithm
}
\description{
Compute rankings based on the transition matrix from a MC algorithm
}
\usage{
MC.ranks(elements, trans, a = 0.15, delta = 10^-15)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{elements}{
Unique elemenets of the union of all input lists - second element of the output list
from function trans.matrix
}
  \item{trans}{
One of the three transition matrices build by function trans.matrix - 4 (5 or 6)
elements of the output list from function trans.matrix
}
  \item{a}{
Tuning parameter to make sure Markov chain with the transition matrix is ergodic;
default set to be 0.15 
}
  \item{delta}{
Convergence criterion for stationary distribution; default set to be 10^-15
}
}
\details{
Compute stationary distribution based on a Markov chain transition matrix built with
function trans.matrix.
}
\value{
Output is a list with 3 components
\item{comp1 }{Number of iteration to reach the stationary distribution}
\item{comp2 }{The stationary distribution}
\item{comp3 }{The rankings based on the stationary distribution}
}
\references{
Lin, S. (2010)
Space oriented rank-based data integration.
Statistical Applications in Genetics and Molecular Biology 9, Article 20.
}
\author{
Shili Lin <shili@stat.osu.edu>
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
trans.matrix
}
\examples{
#get sample data
data(TopKSpaceSampleInput)

mm=trans.matrix(input,space) #platform dependent
MC1=MC.ranks(mm[[2]],mm[[4]])
MC2=MC.ranks(mm[[2]],mm[[5]])
MC3=MC.ranks(mm[[2]],mm[[6]])

mm1=trans.matrix(input,input) #top-k spaces
MC11=MC.ranks(mm1[[2]],mm1[[4]])
MC21=MC.ranks(mm1[[2]],mm1[[5]])
MC31=MC.ranks(mm1[[2]],mm1[[6]])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
