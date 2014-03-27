\name{MC.plot}
\alias{MC.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot of the ordered stationary probabilities
}
\description{
Plot of the ordered stationary probabilities versus ranking contains useful information
regarding the relative rankings of elements. 
}
\usage{

MC.plot(outMC, k, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{outMC}{
Output (a list) from running the MC algorithm 
}
  \item{k}{
Number of stationary probabilities to be plotted. If missing, all stationary probabilities contained in the \code{MC} output will be plotted
}
\item{...}{
Other parameters passed on the the plot function
}
}
\value{
A plot of ordered stationary probabilities versus ranking.
}
\author{
Shili Lin <shili@stat.osu.edu
}

\seealso{
\code{\link{MC}} 
}
\examples{
#get sample data
data(TopKSpaceSampleInput)
outMC=MC(input,space)
MC.plot(outMC)
}
