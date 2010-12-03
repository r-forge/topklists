\name{j0.multi}
\alias{j0.multi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Method for counting K(j0_est)
}
\description{
Method to count K(j0_est) for each pair of several ranked lists. 
}
\usage{
j0.multi(Data,delta, v)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{
input matrix, each coloum presents one list
}
  \item{delta}{
the maximal distance of the objects in ranking
}
  \item{v}{
parameter for estimating j0
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list containing the maximal estimated J0 through all combinations of 2 lists:
\item{maxK}{maximal estimated J0 through all combinations of 2 lists}
\item{K}{matrix with all pair-wise comparison}
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
set.seed(4657)

lists <- data.frame(L1=c("A","B","C","D","E","F","G","H","J","I","K","L","M","N"))
lists$L2 <- c("B","C","A","E","G","F","G","J","K","L","M","N","I","H")
lists$L3 <- sample(LETTERS[1:14])
res.j0.temp = j0.multi(lists, delta=5, v=3)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ list }

