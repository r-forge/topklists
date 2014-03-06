\name{kendallMLists}
\alias{kendallMLists}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  kendallMLists
}
\description{
  Compute Kendall's tau criterion
}
\usage{
kendallMLists(input, aggregate, space = NULL, p = 0.5, w = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{
    A list with each element being a top-\emph{k} list to be aggregated; the top-\emph{k} lists can be of variable lengths
  }
  \item{aggregate}{
    The aggregate list (result) from any of the three classes of algorithms
  }
  \item{space}{
    A list with each element being the underlying space from which the corresponding top-\emph{k} list is derived
  }
  \item{p}{
    A parameter between 0 and 1 for setting the distance of a pair of elements between two lists,  
    if at least one of the elements is not in the underlying space of one
    of the list or if both elements belong to one list but neither belong 
    to the other list. (We recommend using \code{p=0.5} for a "neutral approach".)
  }
  \item{w}{
    Weight vector assigning a weight to each list 
  }
}
\value{
  Kendall's distance
}
\references{
  Lin, S., Ding, J. (2009)
  Integration of ranked lists via Cross Entropy Monte Carlo with applications to mRNA and microRNA studies. Biometrics 65, 9-18.

  Lin, S. (2010)
  Space oriented rank-based data integration. Statistical Applications in Genetics and Molecular Biology 9, Article 20.
}
\author{
  Shili Lin <shili@stat.osu.edu> 
}
\examples{

}

