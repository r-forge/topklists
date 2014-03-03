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
plotBorda(outBorda,K,xlab,ylab,main,sub)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{outBorda}{
A list containing the output from runing the Borda function
}
  \item{K}{

The number of scores to be plotted. If not supplied, all the scores in the output from Borda will be plotted.
}
  \item{xlab}{

Desired text for x-label
}
  \item{ylab}{
Desired text for y-label
}
  \item{main}{
Desired text for title
}
  \item{sub}{
Desired text for sub in plotting function - useful if there are multiple plots in a single figure
}
}
\value{
A plot of Borda's scores versus ranks
}
\references{

}
\author{
 SHili Lin <shili@stat.osu.edu>
}
\seealso{
\code{\link{Borda}}
}
\examples{
#get sample data
data(TopKSpaceSampleInput)
outBorda=Borda(input,space,K=40)
plotBorda(outBorda, K=40)
}

