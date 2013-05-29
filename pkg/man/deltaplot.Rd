\name{deltaplot}
\alias{deltaplot}
\title{
An exploratory plot of discordance for delta selection.
}
\description{
Returns a graph of non-overlap (discordance) of rankings represented by the sum of zeros across all objects in the \eqn{\delta}-dependent Idata vector (see \code{compute.stream}) for a suitable range of \eqn{\delta} values starting at \eqn{\delta=0}. Graphs are plotted for all pairwise list combinations. 
}
\usage{
deltaplot(lists, mind=0, maxd=NULL, subset.lists=NULL, subplot = FALSE, perc.subplot=50) 
}
\arguments{
  \item{lists}{
A data frame cotaining two or more columns that represent lists of ordered objects to be compared
}
  \item{mind}{
Minimal \eqn{\delta} value to be examined; default is 0
}
  \item{maxd}{
Maximal \eqn{\delta} value to be examined, defaults to NULL. If not specified or if \code{maxd} is larger than \code{nrow(lists)}, then \code{maxd=nrow(lists)*0.25}
}
\item{subset.lists}{

}
\item{subplot}{

}
\item{perc.subplot}{
Percentage of the range of the main plot used for creating a subplot in the top right corner, default is 50(\%).  Subplot provides a detailed view of the main plot
}
}
\value{
 \item{Mdelta}{A list of \eqn{\delta}-matrices for each comparison of ordered lists}
}
\author{
 Eva Budinska <budinska@iba.muni.cz>, Michael G. Schimek <michael.schimek@medunigraz.at>, Vendula Svendova <vendula.svendova@medunigraz.at>
}
\references{Schimek, M. G. and Budinska, E. (2010). Visualization techniques for the integration of rank data. In Lechevallier, Y. and Saporta, G. (eds). COMPSTAT 2010. Proceedings in Computational Statistics. Heidelberg: Physica (e-book ISBN 978-3-7908-2603-6), 1637-1644.
}

\examples{
set.seed(1234)
data(TopKGUISampleInput)
##plot subsets
deltaplot(lists, mind=0, maxd=50, subplot=TRUE)

##don't plot subets (default)
deltaplot(lists, mind=0, maxd=50, subplot=FALSE)
}
\keyword{ list}

