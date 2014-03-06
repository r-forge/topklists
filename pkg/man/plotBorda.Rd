\name{plotBorda}
\alias{plotBorda}
\title{
Plot Borda's scores against ranks
}
\description{
 Plotting Borda's scores against ranking can frequently reveal when information for
ranking starts to diminish. This function plots scores versus ranks after aggregation.
}
\usage{
plotBorda(outBorda, k, xlab="Ranking", ylab="Borda Score", main="", sub="")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{outBorda}{
A list containing the output from running the Borda function.
}
  \item{k}{
The number of scores to be plotted. If not supplied, all the scores in the output from Borda will be plotted.
}
  \item{xlab}{
Desired text for x-label.
}
  \item{ylab}{
Desired text for y-label.
}
  \item{main}{
Desired text for title.
}
  \item{sub}{
Desired text for sub in plotting function - useful if there are multiple plots in a single figure.
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
plotBorda(outBorda, k=40)
}

