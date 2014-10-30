\name{deltaplot}
\alias{deltaplot}
\title{
An exploratory plot of discordance for delta selection.
}
\description{
Returns a graph of non-overlap (discordance) of rankings represented by the sum of zeros across all objects in the \eqn{\delta}-dependent \code{Idata} vector (see \code{compute.stream}) for a suitable range of \eqn{\delta} values starting at \eqn{\delta=0}. Graphs are plotted for all pairwise list combinations. 
}
\usage{
deltaplot(lists, deltas=NULL, subset.lists=NULL, subplot = FALSE,
perc.subplot=50, directory=NULL) 
}
\arguments{
  \item{lists}{
A data frame cotaining two or more columns that represent lists of ordered objects to be compared
}
   \item{deltas}{
The range of \eqn{\delta} values to be examined, defaults to NULL. If not specified then \code{delta=c(1:nrow(lists)*0.25)}. If \code{max(deltas)} is larger than \code{nrow(lists)}, then \code{deltas=deltas[which(deltas<nrow(lists)*0.25)]}, 
}
\item{subset.lists}{
Specifies the subset of the input lists, which is used for calculating zero counts for the deltaplot. The value contained in \code{subset.lists} specifies which objects are taken for the calculation from each input list, e.g. a value of 100 would use the first 100 objects of each input list. Default is NULL, in which case all objects of each list are used. 
If specified and \code{max(deltas)} is larger than \code{subset.lists}, then \code{deltas=deltas[which(deltas<subset.lists*0.25)]}, 
}

\item{subplot}{
Logical: if TRUE an additional deltaplot is generated containing a detailed subplot positioned in the top right corner of the plot. This subplot encloses a configurable subset of the values of the original deltaplot. This subset can be specified via a percentage value using the \code{perc.subplot} parameter (default is FALSE for \code{subplot}).
}
\item{perc.subplot}{
Percentage of the range of the main plot used for creating a subplot in the top right corner, default is 50(\%). Subplot provides a detailed view of the main plot.
}
\item{directory}{
Specifies the directory for saving the generated deltaplots in PDF
format. In case \code{directory} is \code{NULL} (default), no pdf is created, but a
new device will be opened.
}


}
\value{
 \item{Mdelta}{A list of \eqn{\delta}-matrices for each comparison of ordered lists}
}
\author{
 Eva Budinska <budinska@iba.muni.cz>, Vendula Svendova <vendula.svendova@medunigraz.at>, Michael G. Schimek <michael.schimek@medunigraz.at>
}
\references{Schimek, M. G. and Budinska, E. (2010). Visualization techniques for the integration of rank data. In Lechevallier, Y. and Saporta, G. (eds). COMPSTAT 2010. Proceedings in Computational Statistics. Heidelberg: Physica (e-book ISBN 978-3-7908-2603-6), 1637-1644.
}

\examples{
set.seed(1234)
data(breast)
##plot subplot
a = deltaplot(breast, deltas = 1:50, subplot=TRUE)

##don't plot subplot (default)
a = deltaplot(breast, deltas=1:50, subplot = FALSE)
}


