\name{prepare.idata}
\alias{prepare.idata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prepare Idata vector of 0's and 1's
}
\description{
Function creates a data stream vector of zeros and ones (\code{Idata}) based on the preselected distance delta of the paired ordered lists. The obtained vector is further used as an input for \code{compute.stream}, a function that estimates the index position of information degradation.
}
\usage{
prepare.idata(x, d)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Data matrix or data frame, where the columns represent the lists of objects ordered according those rankings obtained from two different assessments. 
}
  \item{d}{
The maximal distance between two lists for a ranked object}
}
\details{
The data stream vector is created as follows: if diff(rank1, rank2) of an individual object is less or equal \code{delta}, then 1 is assigned; otherwise 0. The smaller the \code{delta} value, the stronger the assumption of concordance for the paired ranked lists. When \code{delta=0}, the condition returns 1 for an object if and only if its rankings in the two lists are identical (the two objects share the same row).
}
\value{
The result is an object of type \code{Idata}, which is a list containing the data stream vector of zeros and ones, and the information about the applied distance \code{delta}
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
prepare.idata(mm, d=10)

# The breast cancer example
data(breast)
Idata1 = prepare.idata(breast[,c(1,3)], d=10)
# or 
Idata2 = prepare.idata(breast[,c(1,2)], d=10)
# compare to
Idata2 = prepare.idata(breast, d=10)
}


