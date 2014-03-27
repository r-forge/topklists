\name{Kendall.plot}
\alias{Kendall.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot of the Kendall Criterion values
}
\description{
Plot of the Kendall Criterion values of aggregate ranked lists; useful for comparing performances of several algorithms.
}
\usage{
Kendall.plot(input, all.aggregates, space = NULL, algorithm = NULL, p =
0.5, w = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{
A list object containing individual ranked lists.
}
  \item{all.aggregates}{
A list comprising of aggregate top-\emph{k} lists from different algorithms to be compared.
}
  \item{space}{
A list containing the underlying spaces. If not explicitly specified, all lists are treated as coming from a common space defined by the union of all input lists.
}
  \item{algorithm}{
A vector listing the names corresponding to the algorithms used to construct the aggregate ranked lists all.aggregates.
}
  \item{p}{
A parameter between 0 and 1 for setting the distance of a pair of elements between two lists, if at least one of the elements is not in the underlying space of one of the list or if both elements belong to one list but neither belong
to the other list. (We recommend using \code{p=0.5} for a "neutral approach".)
}
  \item{w}{
Weight vector assigned to the input list. Prior information on the reliability of each input list can be incorporated. The default is set to equal weight for each input list.
}
\item{...}{
Other parameters passed on to the plot function.
}

}
\details{
Compute the weighted Kendall's distance between each of the aggregate ranked list with the input ranked lists and plot the computed distances.
}
\value{
A plot of Kendall's distance for each of the aggregate list.
}
\author{
Shili Lin <shili@stat.osu.edu
}

\seealso{
\code{\link{KendallMLists}}
}
\examples{
#get sample data
data(TopKSpaceSampleInput)
outMC=MC(input,space)
all.aggregate=list(outMC$MC1.TopK,outMC$MC2.TopK,outMC$MC3.TopK)
Kendall.plot(input, all.aggregate,space, algorithm=c("MC1","MC2","MC3"))
}
