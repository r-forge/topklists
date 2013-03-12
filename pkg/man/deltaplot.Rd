\name{deltaplot}
\alias{deltaplot}
\title{
An exploratory plot of discordance for delta selection.
}
\description{
Returns a graph of non-overlap (discordance) of rankings represented by the sum of zeros across all objects in the delta-dependent Idata vector (see \code{compute.stream}) for a suitable range of delta values starting at delta=0. Graphs are plotted for all pairwise list combinations. 
}
\usage{
deltaplot(lists, mind=0, maxd=NULL) 
}
\arguments{
  \item{lists}{
a data frame cotaining two or more columns that represent lists of ordered objects to be compared
}
  \item{mind}{
minimal delta value to be examined, default is 0
}
  \item{maxd}{
maximal delta value to be examined, defaults to NULL. If not specified or if maxd is larger than nrow(lists), then maxd=nrow(lists)/2
}
}
\value{
 \item{Mdelta}{a list of delta matrices - for each comparison of ordered lists}
}
\author{
 Eva Budinska <budinska@iba.muni.cz>, Michael G. Schimek <michael.schimek@medunigraz.at>
}
\references{Schimek, M. G. and Budinsk\'{a}, E. (2010). Visualization techniques for the integration of rank data. In Lechevallier, Y. and Saporta, G. (eds). COMPSTAT 2010. Proceedings in Computational Statistics. Heidelberg: Physica (e-book ISBN 978-3-7908-2603-6), 1637-1644.
}

\examples{
set.seed(1234)
data(TopKGUISampleInput)
deltaplot(lists, mind=0, maxd=50)
}
\keyword{ list}

