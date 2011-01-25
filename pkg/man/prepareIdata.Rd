\name{prepareIdata}
\alias{prepareIdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function takes as an input two rankings of a gene list and creates a signal-noise vector of 1's and 0's based on the selected distance delta between rankings. This vector is further used as an input for compute.stream  - a function that estimates the position of signal degeneration
}
\description{
function to prepare Idata from several rankings of several assessors
}
\usage{
prepareIdata(x, delta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Data matrix, where columns represents rank order of genes from two different assessors, rows represents genes (clones)
}
  \item{delta}{
The maximal distance in ranking of the object between two lists }
}
\details{
The signal-noise vector is created as follows: if diff(rank1, rank2) of the object (gene) g is smaller or equal delta, than 1; 0 otherwise. The smaller the delta, the stronger the assumptions on the ressemblence of the two rankings. If delta=0, the condition returns 1 for an object if and only if its rankings in two lists are identical.
}
\value{
The result is a object of type "Idata", which is a list containing the signal-noise vector of 1's and 0's (Idata) and the information
about the distance used (delta).
\item{Idata}{signal-noise vector of 1's and 0's}
\item{delta}{the given delta}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Eva Budinska <budinska@iba.muni.cz>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \link{compute.stream}
}
\examples{
set.seed(4568)
A <- sample(1:20, 20)
B <- sample(1:20, 20)
C <- sample(1:20, 20)
mm <- data.frame(A, B, C, row.names=LETTERS[1:20])
prepareIdata(mm, delta=10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ list }

