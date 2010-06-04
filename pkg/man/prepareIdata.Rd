\name{prepareIdata}
\alias{prepareIdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prepare data from rankings
}
\description{
function to prepare Idata from several rankings of several assessors
}
\usage{
prepareIdata(x, type, genelist, expr = "FALSE", m.dist)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a data matrix containing ranks
}
  \item{type}{
type = first: first column of x represents rank order of genes from the first assessor, second column represents assignment to the category (0, 1; provided by second assessor)
if type = second:  columns of x represents rank order of genes from two different assessors, rows represents genes (clones)
}
  \item{genelist}{
Vector containing clone names. Must be the same length as the number of rows in x.  Order of genes in x must be the same as in genelist. 
}
  \item{expr}{
flag, indicating whether the data are expression values (will be needed in further version of the algorithm) 
}
  \item{m.dist}{
the maximal distance of the objects in ranking between assessors
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
the result is a object of type "Idata", which is a list containing Idata and the information
about the distance (m.dist).
\item{Idata}{ranking ditance information}
\item{m.dist}{distance information}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
set.seed(4568)

A <- sample(1:100, 100)
B <- sample(1:100, 100)
C <- sample(1:100, 100)

prepareIdata(cbind(A, B, C), type="1", genelist=LETTERS[1:10], expr = "FALSE", m.dist=10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ list }

