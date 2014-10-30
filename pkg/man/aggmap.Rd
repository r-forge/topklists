\name{aggmap}
\alias{aggmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Aggregation map for the integration of truncated lists
}
\description{
The function \code{aggmap} applies Paul Murrell's \code{grid} package. It is plotting the ranked items (objects) of \emph{L} truncated (top) lists of individual length \eqn{\hat{k}_i}, based on pairwise comparison of all \emph{L} lists. The resulting \code{aggmap} is defined as follows: For an index, \eqn{p=1,2,\ldots$, $L-1}, aggregation levels (groupings of top lists) are combined in one display. For each group of \emph{L-p} truncated lists down to the smallest group consisting of just one pair of lists, (1) an arbitrary reference list ("ground truth") is selected under the condition that it comprises \eqn{\max_i(\hat{k}_i)} items among all pairwise comparisons in the group of rankings, (2) symbols of its \eqn{\max_i(\hat{k}_i)} items are printed vertically from the highest to the lowest rank position, and (3) the aggregation information for all remaining \emph{L-p} rankings in the group is added, ordered according to descending list length.

The aggregation information per item and group consists of three measures represented by colored triangles and rectangles, respectively, outlined in array format: (a) The \bold{membership} of an individual item in the top-\emph{k} lists, \emph{yes} is denoted by the color 'grey' and \emph{no} by the color 'white'. (b) The \bold{distance} \emph{d} of the rank of an individual item in the reference list from its position in the other list, is denoted by a triangle color scaled from 'red' \emph{identical} to 'yellow' \emph{far distant}. An additional integer value gives the numerical distance between the item's rank positions, a negative sign means ranked lower, and a positive sign means ranked higher in the current list, with respect to the reference list. (c) The rectangular of a symbol takes on the color 'grey' when the \bold{percentage} of \eqn{d \leq \delta} across the columns of a group is above some prespecified threshold, and 'white' otherwise. 
}

\usage{
aggmap(truncated.lists)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{truncated.lists}{
Object resulting from the \code{\link{calculate.maxK}} function
}
}
\author{
Eva Budinska <budinska@iba.muni.cz>, Michael G. Schimek <michael.schimek@medunigraz.at>
}
\references{Murrell, R. (2005) R Graphics. Chapman & Hall/CRC, Boca Raton, Florida.

Schimek, M. G. and Budinska, E. (2010). Visualization techniques for the integration of rank data. In Lechevallier, Y. and Saporta, G. (eds). COMPSTAT 2010. Proceedings in Computational Statistics. Heidelberg: Physica (e-book ISBN 978-3-7908-2603-6), 1637-1644.

Schimek, M. G. and Bloice, M. (2012). Modelling the rank order of Web search engine results. In Komarek, A. and Nagy, S. (eds). Proceedings of the 27th International Workshop on Statistical Modelling. (e-book ISBN 978-80-263-0250-6), Vol. 1, 303-308. 
}
\seealso{
  \code{\link{calculate.maxK}}
}

\examples{
set.seed(1234)
data(breast)
truncated.lists = calculate.maxK(breast, d=10, v=10, L=3, threshold=50)
\dontrun{
aggmap(truncated.lists)
}
}

