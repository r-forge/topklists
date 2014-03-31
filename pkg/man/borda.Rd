\name{Borda}
\alias{Borda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Borda based rank aggregation
}
\description{
Computes Borda scores and ranks based on four different aggregation functions.
}
\usage{
Borda(input, space = NULL, k = NULL)
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
}
\details{
Computes Borda scores and ranks based on four different aggregation functions, in which the underlying spaces, where the individual ranked lists come from, are taken into account. The four aggregation functions are mean, median, geometric mean, and L2 norm.
}
\value{
A list with two components:
\item{TopK}{A matrix with 4 columns each corresponding to the rankings
by each of the 4 aggregation functions.}
\item{Scores}{A matrix with 4 columns each corresponding to the Borda
scores from each of the 4 aggregation functions}
}
\references{
Lin, S. (2010). Space oriented rank-based data integration.
Statistical Applications in Genetics and Molecular Biology 9, Article 20.
}
\author{
Shili Lin <shili@stat.osu.edu>
}
\seealso{
  \link{geo.mean}, \link{l2norm}
  }
\examples{
#get sample data
data(TopKSpaceSampleInput)

outBorda=Borda(input,space) #underlying space-dependent
outBorda1=Borda(input,space=input) #top-k spaces

}

