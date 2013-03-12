\name{j0.multi}
\alias{j0.multi}
\title{
Function returning an overall point j0 of degeneration into noise for multiple ranked lists, combining functions compute.stream and prepare.Idata
}
\description{
Moderate deviation-based calculation of an overall point j0 of degeneration into noise for multiple ranked lists. The function takes a matrix of ordered lists and estimates a j0 for each pair of the lists (columns), with repect to the preselected distance parameter delta.
}
\usage{
j0.multi(lists, d, v)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lists}{
Input matrix, where each column represents one list of ordered items (e.g. gene names ordered according to some measure of effect size) 
}
  \item{d}{
The maximal distance of an object's rank positions when two lists are compared. When the distance between the respective rank positions of the object is smaller or equal d, then the object is assigned the value 1, otherwise 0.
}
  \item{v}{
Parameter for estimating j0
}
}
\details{
The smaller the delta, the stronger assumption on the concordance of any two lists, with delta=0 representing the identity of two rankings of the same object.
}
\value{
A list containing the maximal estimated indices of information degradation j0 through all combinations of L lists:
\item{maxK}{maximal estimated j0 through all combinations of 2 lists}
\item{K}{matrix with estimated j0 for each pair-wise comparison}
}
\author{
Eva Budinska <budinska@iba.muni.cz>, Michael G. Schimek <michael.schimek@medunigraz.at>
}
\references{Hall, P. and Schimek, M. G. (2012). Moderate deviation-based inference for random degeneration in paired rank lists. J. Amer. Statist. Assoc., 107, 661-672.
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
res.j0.temp = j0.multi(lists, d=5, v=3)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ list }

