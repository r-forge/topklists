\name{TopKLists-package}
\alias{TopKLists}
\docType{package}
\title{
Inference, aggregation and visualization for top-\emph{k} ranked lists
}
\description{Web search engines or microarray laboratory devices, among other new technologies, produce very long lists of distinct items or objects in rank order. The statistical task is to identify common top-ranking objects from two or more lists and to form sublists of consolidated items. In each list, the rank position might be due to a measure of strength of evidence, to a preference, or to an assessment either based on expert knowledge or a technical device. For each object, it is assumed that its rank assignment in one list is independent of its rank assignments in the other lists. The ranking is from 1 to \emph{N} throughout without ties. For a general definition of ranked lists see Schimek (2011).

Starting with the work of Mallows (1957), there is a substantial model-based literature on problems in combining rankings where the number of items \emph{N} is relatively small, and significantly less than the number \emph{L} of assessors (rankings). These well-known parametric approaches cannot handle data of the type described above with \emph{N>>L} and \emph{N} huge. Dwork et al. (2001) and DeConde et al. (2006) were the first to address such large-scale rank aggregation problems in the context of Web search engine technology and high-throughput biotechnology, repectively. Here our task is not limited to the aggregation of rankings, we also consider the problem of ranked lists where the reliability of rankings breaks down after the first (top) \emph{k} objects due to error or lack of discriminatory information. In response to the above requirements, we have implemented various distribution-free, and at the same time computationally highly efficient, stochastic approaches because list consolidation by means of brute force (e.g. combinatorial approaches) is limited to the situation where both \emph{N} and \emph{L} are impractically small.

For multiple full ranked (input) lists representing the same set of \emph{N} objects, the package \verb{TopKLists} offers (1) statistical inference on the lengths of informative (top-\emph{k}) partial lists, (2) stochastic aggregation of full or partial lists, and (3) graphical tools for the statistical exploration of input lists, and for aggregation visualization. Our implementations are based on recently developed methods as outlined in Hall and Schimek (2012), Lin (2010a), Lin and Ding (2009), and Schimek, Mysickova and Budinska (2012). Whenever you use the package, please refer to Hall and Schimek (2012) and Lin and Ding (2009) and Schimek et al. (2015) (for full citation please see below).}

\details{
The package consists of three modules and a graphical user interface
(GUI):
\itemize{
\item (1) \verb{TopKInference} provides exploratory nonparametric inference for the estimation of the top-\emph{k} list length of paired rankings;
\item (2) \verb{TopKSpace} provides several rank aggregation techniques (Borda, Markov chain, and Cross Entropy Monte Carlo) which allow the combination of input lists even when the rank positions of some objects are not present in all the lists (so-called \emph{partial input lists});
\item (3) \verb{TopKGraphics} provides a collection of graphical tools for visualization of the inputs to and the outputs from the other modules.
}

Highly convenient is a new aggregation mapping tool called \code{aggmap} in \verb{TopKGraphics}. The GUI allows the non-statistician an easy access to the practically most relevant techniques provided in \verb{TopKInference}, \verb{TopKSpace}, and \verb{TopKGraphics}. Due to the exploratory nature of the implemented methods, tuning parameters are required. All those having a strong impact on the results can be controlled via the GUI. For additional program details and a bioscience application see Schimek et al. (2011). For aspects of modelling the rank order of Web search engine results see Schimek and Bloice (2012). A Springer monograph by Schimek, Lin and Wang of the title ``Statistical Integration of Omics Data'' is in preparation.} 

\author{
Michael G. Schimek, Eva Budinska, Jie Ding, Karl G. Kugler, Vendula Svendova, Shili Lin. 
}
\seealso{
Project homepage: \url{http://topklists.r-forge.r-project.org}
}

\references{
DeConde R. et al. (2006). Combining results of microarray experiments: a rank aggregation approach. Statist. Appl. Genet. Mol. Biol., 5, Article 15.

Dwork, C. et al. (2001). Rank aggregation methods for the Web. \url{http://www10.org/cdrom/papers/577/}

Hall, P. and Schimek, M. G. (2012). Moderate deviation-based inference for random degeneration in paired rank lists. J. Amer. Statist. Assoc., 107, 661-672.

Lin, S. (2010a). Space oriented rank-based data integration. Statist. Appl. Genet. Mol. Biol., 9, Article 20.

Lin, S. (2010b). Rank aggregation methods. Wiley Interdisciplinary Reviews: Computational Statistics, 2, 555-570.

Lin, S. and Ding, J. (2009). Integration of ranked lists via Cross Entropy Monte Carlo with applications to mRNA and microRNA studies. Biometrics, 65, 9-18. 

Mallows, C. L. (1957). Non null ranking models I. Biometrika, 44, 114-130.

Schimek, M. G. (2011). Statistics on Ranked Lists. In Lovric, M. (ed). International Encyclopedia of Statistical Science. Berlin: Springer, Part 19, 1487-1491, DOI: 10.1007/978-3-642-04898-2_563.

Schimek, M. G. and Bloice, M. (2012). Modelling the rank order of Web search engine results. In Komarek, A. and Nagy, S. (eds). Proceedings of the 27th International Workshop on Statistical Modelling. (e-book ISBN 978-80-263-0250-6), Vol. 1, 303-308. 

Schimek, M. G. and Budinska, E. (2010). Visualization Techniques for the Integration of Rank Data. In Lechevallier, Y. and Saporta, G. (eds). COMPSTAT 2010. Proceedings in Computational Statistics. Heidelberg: Physica (e-book ISBN 978-3-7908-2603-6), 1637-1644.

Schimek, M. G., Budinska, E., Kugler, K. and Lin, S. (2011). Package ``TopKLists'' for rank-based genomic data integration. Proceedings of CompBio 2011, 434-440, DOI: 10.2316/P.2011.742-032.

Schimek, M. G., Budinska, E., Kugler, K. G., Svendova, V., Ding, J., Lin, S. (2015). TopKLists: a comprehensive R package for statistical inference, stochastic aggregation, and visualization of multiple omics ranked lists. Statistical Applications in Genetics and Molecular Biology, 14(3): 311-316.

Schimek, M. G., Mysickova, A. and Budinska, E. (2012). An inference and integration approach for the consolidation of ranked lists. Communications in Statistics - Simulation and Computation, 41:7, 1152-1166. 

Schimek, M. G., Lin, S. and Wang, N. (2015). Statistical Integration of Omics Data. In preparation. New York: Springer.
}


