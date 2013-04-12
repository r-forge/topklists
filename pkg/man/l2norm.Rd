\name{l2norm}
\alias{l2norm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the L2 norm
}
\description{
Calculated the L2 norm.
}
\usage{
l2norm(x, na.rm = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
 Objects for which the L2 norm is to be calculated
}
  \item{na.rm}{
Whether or not to remove NA values from the calculation
}
}
\value{
The L2 norm of \code{x}
}
\author{
Shili Lin <shili@stat.ohio-state.edu>
}
\examples{
set.seed(122)
vals <- sample(1:100, 10)
l2norm(vals)
}

