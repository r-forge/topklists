\name{geo.mean}
\alias{geo.mean}
\title{
Calculate the geometric mean
}
\description{
Calculate the geometric mean
}
\usage{
geo.mean(x, na.rm = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A vector of values
}
  \item{na.rm}{
Whether missing values should be automatically removed from calculation
}
}
\value{
The geometric mean
}
\references{
Lin, S. (2010)
Space oriented rank-based data integration.
Statistical Applications in Genetics and Molecular Biology 9, Article 20.
}
\author{
Shili Lin <shili@stat.osu.edu>
}
\seealso{
\link{Borda}
}
\examples{
set.seed(122)
vals <- sample(1:100, 10)
geo.mean(vals)
}

