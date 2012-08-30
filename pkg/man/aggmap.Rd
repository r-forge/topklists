\name{aggmap}
\alias{aggmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Heatmap solution to HS rank order graphics representation
}
\description{
Function plotting the ranked genes in top K lists, based on pairwise comparison of K lists. The resulting graphs is a heatmap of top K ranked genes for each of the pariwise comparisons, starting with the comparison with the longest top K gene list, ordered according to the first list in the comparison. Each cell contains two basic informations: 
 - top right triangle represents information about the distance in ranking of the gene between the two lists represented by number and color code (from red - 0 to white - the maximal distance);
- left bottom triangle represents information on whether the gene from the first list is actually present in the top second list and thus in the top-K ranking (white if yes, grey if not).

}
\usage{
aggmap(lists, res.j0.temp, N, K, d, v)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lists}{
data.frame containing the lists of ranked features, the same as the input to top K ranking functions
}
  \item{res.j0.temp}{
data frame, result of HS algorithm for all pairs of lists, 4 columns, first and second column contain lists names that were compared, third column selected nu and fourth column estimated j0
}
  \item{N}{
Maximal number of features
}
  \item{K}{
Number of lists
}
  \item{d}{
Maximal distance between feature ranks used for j0 estimation
}
  \item{v}{
nu value used for j0 estimation
}
}
\author{
Eva Budinska <budinska@iba.muni.cz>
}
\examples{
set.seed(1234)
lists <- data.frame(L1=c("A","B","C","D","E","F","G","H","J","I","K","L","M","N"))
lists$L2 <- c("B","C","A","E","G","F","D","J","K","L","M","N","I","H")
lists$L3 <- sample(LETTERS[1:14])
res.j0.temp = j0.multi(lists, d=2, v=4)
aggmap(lists, res.j0.temp, N=4, K=3, d=5, v=3)
}
\keyword{ hplot }
