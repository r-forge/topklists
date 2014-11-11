\name{calculate.maxK}
\alias{calculate.maxK}
\title{
The main function for TopKInference
}
\description{
Returns a complex object named truncated.lists containing the \code{Idata} vector (see \code{prepare.idata}), the estimated truncation index \eqn{j_0=k+1} (see \code{compute.stream}) for each pair of input lists, the overall top-\emph{k} estimate (see \code{j0.multi}), and other objects with necessary plotting information for the \code{aggmap}
}
\usage{
calculate.maxK(lists, L, d, v, threshold) 
}
\arguments{
  \item{lists}{
Data frame containing two or more columns that represent input lists of ordered objects subject to comparison
}
  \item{L}{
Number of input lists that are compared
}
  \item{d}{
The maximal distance delta between object ranks required for the estimation of \eqn{j_0}
}
  \item{v}{
The pilot sample size (tuning parameter) \eqn{\nu} required for the estimation of \eqn{j_0}
}
  \item{threshold}{
The percentage of occurencies of an object in the top-\emph{k} selection among all comparisons in order to be gray-shaded in the \code{aggmap} as a consolidated object
}
}
\value{
A named list of the following content:
 \item{comparedLists}{Contains information about the overlap of all pairwise compared lists (structure for the \code{aggmap})}
 \item{info}{Contains information about the list names}
 \item{grayshadedLists}{Contains information which objects in a list are consolidated (gray-shaded in the \code{aggmap})}
 \item{summarytable}{Table of top-\emph{k} list overlaps containing rank information, the rank sum, the order of objects as a function of the rank sum, the frequency of an object in the input lists and the frequency of an object in the truncated lists (for plotting in the \code{aggmap})}
 \item{vennlists}{Contains the top-\emph{k} objects for each of the input lists (for display in the Venn-diagram)}
 \item{venntable}{Contains the overlap information (for display in the Venn-table)}
 \item{v}{Selected pilot sample size (tuning parameter) \eqn{\nu}}
 \item{Ntoplot}{Number of columns to be plotted in the \code{aggmap}}
 \item{Idata}{Data frame of Idata vectors (see \code{compute.stream}) for each pair of input lists and the associated delta's}
 \item{d}{selected delta}
 \item{threshold}{selected threshold}
 \item{threshold}{number of lists}
 \item{N}{number of items in data frame (lists)}
 \item{lists}{data frame of lists that entered the analysis}
 \item{maxK}{maximal estimate of the top-\emph{k}'s (for all pairwise comparisons)}
 \item{topkspace}{the final integrated list of objects as result of the CEMC algorithm applied to the maxK truncated lists}
 
}
\author{
 Eva Budinska <budinska@iba.muni.cz>, Michael G. Schimek <michael.schimek@medunigraz.at>
}

\references{
Hall, P. and Schimek, M. G. (2012). Moderate deviation-based inference for random degeneration in paired rank lists. J. Amer. Statist. Assoc., 107, 661-672.
}
\seealso{
  \code{\link{CEMC}, \link{prepare.idata}}
}
\examples{
set.seed(1234)
data(breast)
truncated.lists = calculate.maxK(breast, d=6, v=10, L=3, threshold=50)
\dontrun{
aggmap(truncated.lists)
}
}


