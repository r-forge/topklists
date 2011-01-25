\name{j0.multi}
\alias{j0.multi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function returning the estimated maximal point of signal degeneration j0 between two or multiple lists, combines functions compute.stream and prepare.Idata
}
\description{
Method to estimate j0 for each pair of several ranked lists. The function takes a matrix of ordered lists and estimates the point of degeneration of signal to noise for each pair of the lists (columns), based on selected maximal distance between two items in lists delta.
}
\usage{
j0.multi(Data,delta, v)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{
input matrix, where each column represents one list of ordered items (e.g. gene names ordered according to increasing p-value) 
}
  \item{delta}{
the maximal distance in rankings of the object between two lists, used for constructing the signal-noise vector (if the distance between two rankings of the object in two lists is smaller or equal delta, than 1, otherwise 0).
}
  \item{v}{
parameter for estimating j0
}
}
\details{
The smaller the delta, the stronger assumption on the ressemblance between any two lists, with delta=0 representing the identity of two rankings of the same object.
}
\value{
A list containing the maximal estimated point of signal degeneration j0 through all combinations of 2 lists:
\item{maxK}{maximal estimated j0 through all combinations of 2 lists}
\item{K}{matrix with estimated j0 for each pair-wise comparison}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Eva Budinska <eva.budinska@isb-sib.ch>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link{compute.stream}, \link{prepareIdata}
}
\examples{
set.seed(4657)

lists <- data.frame(L1=c("A","B","C","D","E","F","G","H","J","I","K","L","M","N"))
lists$L2 <- c("B","C","A","E","G","F","G","J","K","L","M","N","I","H")
lists$L3 <- sample(LETTERS[1:14])
res.j0.temp = j0.multi(lists, delta=5, v=3)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ list }

