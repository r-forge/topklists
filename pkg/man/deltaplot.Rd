\name{deltaplot}
\alias{deltaplot}
\title{
An exploratory plot of discordance for delta selection.
}
\description{
Returns a graph of non-overlap (discordance) of rankings represented by the sum of zeros across all objects in the \eqn{\delta}-dependent Idata vector (see \code{compute.stream}) for a suitable range of \eqn{\delta} values starting at \eqn{\delta=0}. Graphs are plotted for all pairwise list combinations. 
}
\usage{
deltaplot(lists, mind=0, maxd=NULL, perc.subplot=50, subset.plotted=NA) 
}
\arguments{
  \item{lists}{
A data frame cotaining two or more columns that represent lists of ordered objects to be compared
}
  \item{mind}{
Minimal \eqn{\delta} value to be examined; default is 0
}
  \item{maxd}{
Maximal \eqn{\delta} value to be examined, defaults to NULL. If not specified or if \code{maxd} is larger than \code{nrow(lists)}, then \code{maxd=nrow(lists)/2}
}
  \item{perc.subplot}{
Maximal \eqn{\delta} value to be examined, defaults to NULL. If not specified or if \code{maxd} is larger than \code{nrow(lists)}, then \code{maxd=nrow(lists)/2}
}
  \item{subset.plotted}{
The default \code{NA} will suppress the creation of any pair-wise subplots. Setting it to \code{NULL} will create all possible pairs of list combinations for creating subplots.
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
##don't make pair-wise comparisons (default)
deltaplot(lists, mind=0, maxd=50, subset.plotted=NULL)

##make all pair-wise comparisons
deltaplot(lists, mind=0, maxd=50, subset.plotted=NULL)
}
\keyword{ list}

