\name{MC}
\alias{MC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Markov chain based rank aggregation
}
\description{
Aggregating ranked lists using three  Markov chain algorithms.
}
\usage{
MC(input, space = NULL, k = NULL, a = 0.15, delta = 10^-15)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{
A list containing individual ranked lists.
}
  \item{space}{
A list containing the underlying spaces. If not explicitly specified, all lists are treated as coming from a common space defined by the union of all input lists.
}
  \item{k}{
An integer specifying the number of items in the output top-k list.
}
  \item{a}{
Tuning parameter to make sure Markov Chain with the transition matrix is ergodic;
default set to 0.15.
}
  \item{delta}{
Convergence criterion for stationary distribution; default set to 10^-15.
}
}
\details{
Constructs ergodic Markov Chain based on ranking data from individual lists. A larger probability in the stationary distribution corresponds to a higher rank of the corresponding element. The algorithm are considered: MC1 (spam sensitive), MC2 (majority rule), and MC3 (proportional).
}
\value{
A list of elements, two for each of the MC algorithms: 
\item{MC1.TopK}{A vector of aggregate ranked elements based on \code{MC1} algorithm.}
\item{MC1.Prob}{Stationary probability distribution: a vector of probabilities corresponding to the ranked elements in \code{MC1.TopK}}
\item{MC2.TopK}{A vector of aggregate ranked elements based on MC2 algorithm.}
\item{MC2.Prob}{Stationary probability distribution: a vector of probabilities corresponding to the ranked elements in \code{MC2.TopK}}
\item{MC3.TopK}{A vector of aggregate ranked elements based on MC3 algorithm.}
\item{MC3.Prob}{Stationary probability distribution: a vector of probabilities corresponding to the ranked elements in \code{MC3.TopK}}
}
\references{
Lin, S. (2010). Space oriented rank-based data integration.
Statistical Applications in Genetics and Molecular Biology 9, Article 20.
}
\author{
Shili Lin <shili@stat.osu.edu>
}
\seealso{
\code{\link{Borda}}, \code{\link{CEMC}}, \code{\link{MC.plot}}
}
\examples{
#get sample data
data(TopKSpaceSampleInput)
outMC=MC(input,space) #underlying space-dependent
outMCa=MC(input,space=input) #top-k spaces
MC.plot(outMC)
}
