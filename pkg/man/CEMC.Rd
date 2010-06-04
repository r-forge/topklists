\name{CEMC}
\alias{CEMC}
\title{
Perform cross-entropy Monte Carlo simulations
}
\description{
Perform cross-entropy Monte Carlo simluations
}
\usage{
CEMC(topK, k, dm = "s", N = NULL, N1 = NULL, rho = 0.1, e1 = 0.1, e2 = 1, w = 0.5, b = 0, init.m = "p", init.w = 0, d.w = NULL, input.par = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{topK}{
A list of several top K lists, may have different length, one column for each list, ended with NA's
}
  \item{k}{
Length of needed list
}
  \item{dm}{
Distance measure, "s" for spearman, "k" for kendall(p=0)
}
  \item{N}{
Number of samples each iterate
}
  \item{N1}{
Number of best samples retained each iterate
}
  \item{rho}{
Proportion of samples used to estimated new probability matrix (default = 0.1)
}
  \item{e1}{
for p
}
  \item{e2}{
for y
}
  \item{w}{
Weight of the new probabilities estimate each iterate
}
  \item{b}{
Parameter used in .blur function
}
  \item{init.m}{
Initialization method, see the function \code{\link{init.p}} for details
}
  \item{init.w}{
Probability matrix initialization:
}
  \item{d.w}{
Weights for distances from different input lists
}
  \item{input.par}{
Dataframe of input parameters
}
}
\details{
}
\value{
CEMC returns names list containing the results, as well as information on the method call
}
\references{
}
\author{
Shili Lin <shili@stat.ohio-state.edu>
}
\note{
}

\seealso{
}
\examples{

A <- sample(1:10, 10)
B <- sample(1:10, 10)
C <- sample(1:10, 10)

CEMC(cbind(A, B, C), 10)
}
\keyword{ robust }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
