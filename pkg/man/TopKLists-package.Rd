\name{TopKLists-package}
\alias{topklists-package}
\alias{TopKLists}
\alias{topklists}
\docType{package}
\title{
Inference, aggregation and visualization for top-k ranked lists
}
\description{Web search engines or microarray laboratory devices, among other new technologies, produce very long lists of distinct items or objects in rank order. The statistical task is to identify common top-ranking objects from two or more lists and to form sublists of consolidated items. In each list, the rank position might be due to a measure of strength of evidence, to a preference, or to an assessment either based on expert knowledge or a technical device. For each object, it is assumed that its rank assignment in one list is independent of its rank assignments in the other lists. The ranking is from 1 to \emph{N} throughout without ties. 

Starting with the work of Mallows (1957), there is a substantial model-based literature on problems in combining rankings where the number of items \emph{N} is relatively small, and significantly less than the number \emph{L} of assessment mechanisms. These well-known parametric approaches cannot handle data of the type described above with \emph{N>>L} and \emph{N} huge. Moreover, here we have problems in mind where the reliability of rankings breaks down after the first (top) \emph{k} objects due to error or lack of discriminatory information. In response to these requirements, we have implemented distribution-free, and at the same time computationally highly efficient, stochastic approaches because list aggregation by means of brute force is limited to the situation where both \emph{N} and \emph{L} are impractically small.

For multiple full ranked (input) lists representing the same set of \emph{N} objects, the package TopKLists offers (1) statistical inference on the lengths of informative (top-\emph{k}) partial lists, (2) stochastic aggregation of full or partial lists, and (3) graphical tools for the statistical exploration of input lists, and for aggregation visualization.}

\details{
The package consists of three modules:
(1) \verb{TopKInference} provides exploratory nonparametric inference for the estimation of the top-\emph{k} list length of paired rankings;
(2) \verb{TopKSpace} provides several rank aggregation techniques (Borda, Markov chain, and Cross-entropy Monte Carlo) which allow the combination of input lists of different lengths (space) as well as missing rank information;
(3) \verb{TopKGraphics} provides a collection of graphical tools for visualization of the inputs to respectively the outputs from the other modules.

For more details and a bioscience application see Schimek et al. (2011).}

\author{
Michael G. Schimek, Eva Budinska, Jie Ding, Peter Hall, Karl G. Kugler, Shili Lin
Maintainers: Michael G. Schimek <michael.schimek@medunigraz.at> and Karl G. Kugler <karl@eigenlab.net>
}
\seealso{
project homepage: \url{http://topklists.r-forge.r-project.org}
}
\references{
DeConde R. P., Hawley, S., Falcon, S., Clegg, N., Knudsen, B. and Etzioni, R. (2006). Combining results of microarray experiments: a rank aggregation approach. Statist. Appl. Genet. Mol. Biol., 5, Article 15.

Hall, P. and Schimek, M. G. (2010). Moderate deviation-based inference for random degeneration in paired rank lists. J. Amer. Statist. Assoc., under revision.

Lin, S. and Ding, J. (2009). Integration of ranked lists via Cross Entropy Monte Carlo with applications to mRNA and microRNA studies. Biometrics, 65, 9-18. 

Lin, S. (2010). Space oriented rank-based data integration. Statist. Appl. Genet. Mol. Biol., 20, Article 20.

Lin, S. (2010). Rank aggregation methods. Wiley Interdisciplinary Reviews: Computational Statistics, 2, 555-570.

Mallows, C. L. (1957). Non null ranking models I. Biometrika, 44, 114-130.

Schimek, M. G. and Budinsk\'{a}, E. (2010). Visualization Techniques for the Integration of Rank Data. In Lechevallier, Y. and Saporta, G. (eds) COMPSTAT 2010. Proceedings in Computational Statistics. Heidelberg: Physica (e-book ISBN 978-3-7908-2603-6), 1637-1644.

Schimek, M. G. (2011). Statistics on Ranked Lists. In Lovric, M. (ed) International Encyclopedia of Statistical Science. Berlin: Springer, Part 19, 1487-1491, DOI: 10.1007/978-3-642-04898-2_563.

Schimek, M. G., Budinsk\'{a}, E., Kugler, K. and Lin, S. (2011). Package ``TopKLists'' for rank-based genomic data integration. Proceedings of the IASTED International Conference Computational Bioscience (CompBio 2011), 434-440, DOI: 10.2316/P.2011.742-032.
}
\keyword{package}

