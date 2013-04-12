\name{Borda}
\alias{Borda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute Borda scores and ranks
}
\description{
Computes Borda scores and ranks based on four different aggregation functions
}
\usage{
Borda(input, space)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{
A list containing individual ranked lists
}
  \item{space}{
A list containing the underlying spaces
}
}
\details{
Computes Borda scores and ranks based on four different aggregation functions, in which the underlying spaces where the individual ranked lists come from, are taken into account. The four aggregation functions are mean, median, geometric mean, and l2 norm.
}
\value{
Output is a list with two components
\item{Comp 1 (RankedObjects)}{A matrix with 4 columns each corresponding to the rankings
by each of the 4 aggregation functions}
\item{Comp 2 (BordaScores) }{A matix with 4 columns each corresponding to the Borda
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
geomean, l2norm}
\examples{
#get sample data
data(TopKSpaceSampleInput)

bb1=Borda(input,space) #platform dependent
bb2=Borda(input,input) #top-k spaces

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }

