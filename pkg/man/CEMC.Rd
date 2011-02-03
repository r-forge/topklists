\name{CEMC}
\alias{CEMC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
CEMC combination of TopKLists
}
\description{
Perform cross-entropy Monte Carlo simluations for generating combined ranked list using CEMC
taking into account difference spaces of input ranked lists
}
\usage{
CEMC(input, k, space = NULL, dm = "Kendall", kp = 0.5, N = NULL, N1 = NULL, rho = 0.1, e1 = 0.1, e2 = 1, w = 0.5, b = 0, init.m = "p", init.w = 0, d.w = NULL, input.par = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{
A list of several TopKLists, may have different length
}
  \item{k}{
Desired length of combined list
}
  \item{space}{
Underlying spaces for the lists, NULL means all lists share common space
}
  \item{dm}{
Distance measure, "Kendall" or "Spearman"
}
  \item{kp}{
Partial distance used in Kendall's tau distance measure
}
  \item{N}{
Number of samples generated in each iterate
}
  \item{N1}{
Number of samples retained after each iterate 
}
  \item{rho}{
Proportion of samples used to estimate new probability matrix
}
  \item{e1}{
Stop criteria w.r.t. p
}
  \item{e2}{
Stop criteria w.r.t. y
}
  \item{w}{
Weight of the new probabilities estimate each iterate
}
  \item{b}{
Parameter used in blur function
}
  \item{init.m}{
Initialization method, see the function \code{init.p} for details
}
  \item{init.w}{
Probability matrix initialization. (See Details)
}
  \item{d.w}{
Weights for distances from different input lists
}
  \item{input.par}{
Input parameters in a data.frame
}
}
\details{
init.w: probability matrix initialization:
(1-init.w) * uniform +  init.w * estimated from input lists
}
\value{
}
\references{
}
\author{
Jie Ding <ding@stat.osu.edu>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
set.seed(1234)
a <- LETTERS[sample(1:20, 10)]
b <- LETTERS[sample(1:20, 10)]
c <- LETTERS[sample(1:20, 10)]

myin <- list(a, b, c)

CEMC(myin, 10, init.m= "p")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
