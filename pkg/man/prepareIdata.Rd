\name{prepareIdata}
\alias{prepareIdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function that reads in two different rankings of a set of objects and creates a data stream vector of zeros and ones based on the preselected distance delta of the paired ranked lists 
}
\description{
Function preparing Idata from several rankings (of several assessors). The obtained vector is further used as an input for \code{compute.stream}, a function that estimates the index position of information degradation.
}
\usage{
prepareIdata(x, d)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Data matrix, where the columns represent the rank order of the objects based on two different assessments, and the rows represent all ranked objects
}
  \item{d}{
The maximal distance between two lists for a ranked object}
}
\details{
The data stream vector is created as follows: if diff(rank1, rank2) of an individual object is less or equal \code{delta}, then 1 is assigned; otherwise 0. The smaller the \code{delta} value, the stronger the assumption of concordance for the paired ranked lists. When \code{delta=0}, the condition returns 1 for an object if and only if its rankings in the two lists are identical.
}
\value{
The result is a object of type "Idata", which is a list containing the data stream vector of zeros and ones (Idata) and the information about the applied distance \code{delta}
\item{Idata}{Data stream vector of zeros and ones}
\item{delta}{The applied \code{delta}}
}
\author{
Eva Budinska <budinska@iba.muni.cz>, Michael G. Schimek <michael.schimek@medunigraz.at>
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
prepareIdata(mm, d=10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ list }

