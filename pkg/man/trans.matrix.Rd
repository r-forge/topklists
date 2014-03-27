\name{trans.matrix}
\alias{trans.matrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute transition matrices
}
\description{
Builds transition matrices for all three Markov Chain algorithms
}
\usage{
trans.matrix(input, space)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{
A list containing individual ranked lists
}
  \item{space}{
A list containing the underlying spaces
}
}
\details{
Both input and space are lists of the same \code{length = nList}
}
\value{
The output is a list: 

\item{L }{Unique elements of the union of all input ranked lists}


\item{MC1 }{The transition matrix constructed from the MC1 algorithm}
\item{MC2 }{The transition matrix constructed from the MC2 algorithm}
\item{MC3 }{The transition matrix constructed from the MC3 algorithm}
}
\references{
Lin, S. (2010)
Space oriented rank-based data integration.
Statistical Applications in Genetics and Molecular Biology 9, Article 20.
}
\author{
Shili Lin <shili@stat.osu.edu>
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{MC}}
}


