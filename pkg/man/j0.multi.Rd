\name{j0.multi}
\alias{j0.multi}
\title{
Function returning an overall point j0 of degeneration into noise for multiple ranked lists
}
\description{
Moderate deviation-based calculation of an overall point \eqn{j_0} of degeneration into noise for multiple ranked lists. The function takes a matrix of ordered lists and estimates a \eqn{j_0} for each pair of the input lists (columns), with repect to the preselected distance parameter \eqn{\delta}. This function combines the functions \code{compute.stream} and \code{prepare.Idata}.
}
\usage{
j0.multi(lists, d, v)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lists}{
Input data frame, where each column represents one list of ordered items
}
  \item{d}{
The maximal distance of an object's rank positions when two lists are compared. When the distance between the respective rank positions of the object is smaller or equal \code{d}, then the object is assigned the value 1, otherwise 0
}
  \item{v}{
Parameter for estimating \eqn{j_0}
}
}
\details{
The smaller \code{d}, the stronger the assumption about the concordance of any two lists (\code{d=0} is assuming identical rankings of an object)
}
\value{
A list containing the maximal estimated indices of information degradation \eqn{j_0} through all combinations of \emph{L} lists:
\item{maxK}{Maximal estimated k through all combinations of two lists}
\item{L}{Data frame of estimated \eqn{j_0} for each pairwise comparison}
\item{Idata}{Data stream vector of zeros and ones}
}
\author{
Eva Budinska <budinska@iba.muni.cz>, Michael G. Schimek <michael.schimek@medunigraz.at>
}
\references{Hall, P. and Schimek, M. G. (2012). Moderate deviation-based inference for random degeneration in paired rank lists. J. Amer. Statist. Assoc., 107, 661-672.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link{compute.stream}, \link{prepare.idata}
}
\examples{
set.seed(4657)

lists <- data.frame(L1=c("A","B","C","D","E","F","G","H","J","I","K","L","M","N"))
lists$L2 <- c("B","C","A","E","G","F","G","J","K","L","M","N","I","H")
lists$L3 <- sample(LETTERS[1:14])
res.j0.temp = j0.multi(lists, d=5, v=3)


}


