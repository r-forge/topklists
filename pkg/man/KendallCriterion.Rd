\name{KendallCriterion}
\alias{KendallCriterion}
\title{
Compute Kendall's criterion
}
\description{
Compute Kendall's criterion
}
\usage{
KendallCriterion(input, space, aggregate, p, w)
}
\arguments{
  \item{input}{
A list with each element being a top-K list to be aggregated. The top-K lists can be of variable lengths.
}
  \item{space}{
A list with each element being the underlying space from which the corresponding top-K list is derived.
}
  \item{aggregate}{
The aggregate list (result) from any of the three classes of algorithms.}
  \item{p}{
A parameter between 0-1 for setting the distance of a pair of elements between two lists  
if at least one of the elements is not in the underlying space of one
of the list or if both elements belong to one list but neither belong 
to the other list. We recommend using p=0.5 for a "neutral approach".
}
  \item{w}{
Weight vector assigning a weight to each list. 
}
}
\value{
A real number containing the Kendall's distance value.
}
\references{
Lin, S., Ding, J. (2009)
Integration of Ranked Lists  via
Cross Entropy Monte Carlo with Applications to mRNA and microRNA
Studies. Biometrics 65, 9-18. <p>
Lin, S. (2010)
Space oriented rank-based data integration.
Statistical Applications in Genetics and Molecular Biology 9, Article 20.
}
\author{
Shili Lin <shili@stat.osu.edu> 
}
\examples{
L1=1:30
L2=c(1:10,31:40,11:15)
L3=c(1:10,16:20,11:15)
input=list(L1,L2,L3)
space=list(1:40,1:40,1:40)
bb1=Borda(input,space) 
kc.ARM=KendallCriterion(input, space, bb1[[1]][, 1], 
p = 0.5, w = c(2/(30 * (30 - 1)), 2/(25 * (25 - 1)), 2/(20 * (20 -+ 1))))
}
\keyword{ univar }
