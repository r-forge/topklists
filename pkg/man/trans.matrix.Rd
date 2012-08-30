\name{trans.matrix}
\alias{trans.matrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute transition matrices
}
\description{
Build transition matrices for all three Markov chain algorithms
}
\usage{
trans.matrix(input, space)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{
a list containing individual ranked lists
}
  \item{space}{
a list containing the underlying spaces
}
}
\details{
Both input and space are lists of the same length=nList
}
\value{
The output is a list 
\item{comp1 }{a list with three components - each component is an extended top-k list
of the same length as the corresponding space}
\item{comp2 }{unique elements of the union of all input ranked lists}
\item{comp3 }{The number of unique element, i.e., length of component 2}
\item{Comp4-6 }{Each of the three components gives the transition matrix constructed
from one of the three Markov chain algorithms
}
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
MC.ranks
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
