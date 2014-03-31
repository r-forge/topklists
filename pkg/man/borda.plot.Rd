\name{Borda.plot}
\alias{Borda.plot}
\title{
Plot Borda's scores against ranks
}
\description{
 Plotting Borda's scores against ranking can frequently reveal when information for
ranking starts to diminish. This function plots scores versus ranks after aggregation.
}
\usage{
Borda.plot(outBorda, k, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{outBorda}{
A list containing the output from running the Borda function.
}
  \item{k}{
The number of scores to be plotted. If not supplied, all the scores in the output from Borda will be plotted.
}
\item{...}{
  other parameters passed on to the plot function
}
}
\value{
A plot of Borda's scores versus ranks.
}
\author{
 Shili Lin <shili@stat.osu.edu>
}
\seealso{
\code{\link{Borda}}
}
\examples{
#get sample data
data(TopKSpaceSampleInput)
outBorda=Borda(input,space,k=40)
Borda.plot(outBorda, k=40)
}

