\name{Borda}
\alias{Borda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute Borda scores and ranks
}
\description{
Compute Borda scores and ranks based on four aggregation function
}
\usage{
Borda(input, space)
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
}
\value{
\item{RankedObjects}{}
\item{BordaScores }{}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Shili Lin <shili@stat.ohio-state.edu>
}
\note{
%%  ~~further notes~~
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#get sample data
data(TopKSpaceSampleInput)

bb1=Borda(input,space) #platform dependent
bb2=Borda(input,input) #top-k spaces

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }

