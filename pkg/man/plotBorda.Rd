\name{plotBorda}
\alias{plotBorda}
\title{
Function to plot Borda's scores and ranks
}
\description{
 Plotting Borda's scores against ranking can frequently reveal when information for
ranking starts to diminish. This function function plots scores and ranks after aggregation.
}
\usage{
plotBorda(input, space, brange)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{
A list containing individual ranked lists
}
  \item{space}{
A list containing the underlying spaces
}
  \item{brange}{
The range to be plotted. Set to the max length of \code{input} lists by default.
}
}
\value{
A plot of Borda's scores and ranks
}
\references{
TODO
}
\author{
 SHili Lin
}
\seealso{
\code{\link{Borda}}
}
\examples{
#get sample data
data(TopKSpaceSampleInput)
plotBorda(input, space)
}

