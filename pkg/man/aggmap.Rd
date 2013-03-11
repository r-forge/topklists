\name{aggmap}
\alias{aggmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Aggregation map for the integration of truncated lists
}
\description{
Function plotting the ranked items of \emph{L} truncated (top) lists of individual length \eqn{\hat{k}_i}, based on pairwise comparison of all \emph{L} lists. The resulting \verb{aggmap} is defined as follows: For an index \eqn{p=1,2,\ldots$, $L-1} aggregation levels (groupings of top lists) are combined in one display. For each group of \emph{L-p} truncated lists down to the smallest group consisting of just one pair of lists, (1) an arbitrary reference list is selected under the condition that it comprises \eqn{\max_i(\hat{k}_i)} items among all pairwise comparisons in the group of rankings, (2) symbols of its \eqn{\max_i(\hat{k}_i)} items are printed vertically from the highest to the lowest rank position, and (3) the aggregation information for all remaining \emph{L-p} rankings in the group is added, ordered according to descending list length.

The aggregation information per symbol, item, and group consists of three measures represented by colored triangles and rectangles, respectively, outlined in array format: (a) The \emph{membership} of an individual item in the top-\emph{k} lists, \emph{yes} is denoted by the color 'grey' and \emph{no} by the color 'white'. (b) The \emph{distance} \emph{d} of the rank of an individual item in the reference list from its position in the other list, is denoted by a triangle color scaled from 'red' \emph{identical} to 'yellow' \emph{far distant}. An additional integer value gives the numerical distance between the item's rank positions, a negative sign means ranked lower, and a positive sign means ranked higher in the current list with respect to the reference list. (c) The rectangular of a symbol takes on the color 'grey' when the \emph{percentage} of \eqn{d \leq \delta} across the columns of a group is above some prespecified threshold, and 'white' otherwise. 
}

\usage{
aggmap(lists, res.j0.temp, N, L, d, v)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lists}{
data.frame containing the lists of ranked objects, the same as the input to top k ranking functions
}
  \item{res.j0.temp}{
data frame, result of truncation algorithm for all pairs of lists, 4 columns, first and second column contain lists names that were compared, third column selected nu and fourth column estimated j0
}
  \item{N}{
Maximal number of objects
}
  \item{L}{
Number of ranked lists
}
  \item{d}{
Maximal distance between object ranks used for j0 estimation
}
  \item{v}{
nu value used for j0 estimation
}
}
\author{
Eva Budinska <budinska@iba.muni.cz>, Michael G. Schimek <michael.schimek@medunigraz.at>
}
\examples{
set.seed(1234)
data(TopKGUISampleInput)
res.j0.temp = j0.multi(lists, d=2, v=4)
\dontrun{
aggmap(lists, res.j0.temp, N=4, L=3, d=5, v=3)
}
}
\keyword{ hplot }
