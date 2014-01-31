\name{plotMC}
\alias{plotMC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot of the ordered stationary probabilities
}
\description{
Plot of the ordered stationary probabilities versus ranking contains useful information
regarding the relative rankings of elements. 
}
\usage{
plotMC(input, space, prob.range)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{
A list containing individual ranked lists
}
  \item{space}{
A list containing the underlying spaces
}
  \item{prob.range}{
The range to be plotted. Set to the max length of \code{input} lists by default.
}
}
\details{
TODO
}
\value{
A plot of ordered stationary probabilities versus ranking.
}
\references{
TODO
}
\author{
Shili Lin
}

\seealso{
\code{\link{trans.matrix}}, \code{\link{MC.ranks}} 
}
\examples{
#get sample data
data(TopKSpaceSampleInput)
plotMC(input, space)
}

