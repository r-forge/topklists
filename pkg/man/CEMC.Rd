\name{CEMC}
\alias{CEMC}
\title{
CEMC based rank aggregation
}
\description{
Performs Cross Entropy Monte Carlo simulations for generating combined ranked list using CEMC, taking into account the different spaces of ranked input lists.
}
\usage{
CEMC(input, space = NULL, k=NULL, dm = "k", kp = 0.5, N = NULL, N1 = NULL,
rho = 0.1, e1 = 0.1, e2 = 1, w = 0.5, b = 0, init.m = "p", init.w = 0,
d.w = NULL, input.par = NULL, extra=0)
}
\arguments{
  \item{input}{
A list of several \code{TopKLists}, may have different length
}
  \item{space}{
A list of the same structure as the input list. Contains underlying spaces for the top-\emph{k} lists. NULL means all lists share a common space, which is taken to be the union of all input lists
}
  \item{k}{
Desired length of combined list
}
  \item{dm}{
Distance measure, "s" for Spearman, "k" for Kendall (p=0)
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
Proportion of samples used to estimate a new probability matrix
}
  \item{e1}{
Stopping criterion with respect to the l1-norm of the difference of the two probability matrices between the current and previous iterations  
}
  \item{e2}{
Stopping criterion with respect to the difference in the obtimizing criterion (e.g. the generalized Kemeny guideline) between the current and the previous iterations  
}
  \item{w}{
Weight of the new probability vector for the next iterate
}
  \item{b}{
Parameter used in blur function - this is for finding starting values for the algorithem
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
  \item{extra}{
Number of additional items to be included in the combined ranked list during the 
calculation
}
}
\details{
The algorithm implemented is the Order Explicit Algorithm, which is an iterative procedure to maximize an objective function (either based on Kendall's distance (dm="k") or Spearman's distance (dm="s")). 

init.w: probability matrix initialization:
(1-init.w) * uniform +  init.w * estimated from input lists
}
\value{
A list containing three components:
\item{TopK}{A vector giving the aggregate ranked list.}
\item{ProbMatrix}{A matrix, with each column represent the probability vector of a multinomial distribution and thus sum to 1.}
\item{input.par}{A vector containing tuning parameters used in the current run. User may edit this vector and use it as input for a more refined analysis.
}
}
\references{
Lin, S. and Ding, J. (2009). Integration of ranked lists via Cross Entropy Monte Carlo with applications to mRNA and microRNA studies. Biometrics, 65, 9-18.
}
\author{
Jie Ding <jding@jimmy.harvard.edu>, Shili Lin <shili@stat.osu.edu>
}
\examples{
#small data set; a larger data example is available in the vignettes
L1=c("chicken","dog","cat")
L2=c(1,"chicken","cat", 2:5)
L3=c("dog","chicken",1:10)
input=list(L1,L2,L3)
space1=c("chicken","dog","cat",1:10)
space=list(space1,space1,space1)
outCEMC=CEMC(input, space) #underlying space-dependent
}
