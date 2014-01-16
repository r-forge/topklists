\name{breast}
\docType{data}
\alias{breast}
\title{Sample data from breast cancer expression}
\description{
The example dataset comprises three lists of microarray results (differential gene
expression) from three breast cancer studies:
\enumerate{
 \item MicroArray Quality Control Phase II Project (2010) labeled MDACC
 \item Strong Time Dependence of the 76-Gene Prognostic Signature (2007) labeled TransBig
 \item A Clinically Relevant Gene Signature in Triple-Negative and Basal-Like Breast Cancer (2011) labeled Pusztai
}
Only genes (unique gene symbols) common to all studies are considered, therefore
each of the three lists has the length of N = 917.
}
\usage{breast}
\format{data.frame}
\source{
\enumerate{
 \item MDACC: \href{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE20194}{GEO GSE20194}
 \item TransBig: \href{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse7390}{GEO GSE7390}
 \item Pusztai: \href{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE31519}{GEP GSE31519}
}
}
\references{
\enumerate{
	\item Desmedt C, Piette F, Loi S, Wang Y et al. (2007).  Strong time dependence of the 76-gene prognostic signature for node-negative breast cancer patients in the TRANSBIG multicenter independent validation series. Clin Cancer Res, 1;13(11):3207-14. PMID: 17545524
	\item Karn T, Pusztai L, Holtrich U, Iwamoto T et al. (2011). Homogeneous datasets of triple negative breast cancers enable the identification of novel prognostic and predictive signatures. PLoS One ;6(12):e28403. PMID: 22220191
	\item Karn T, Pusztai L, Ruckhäberle E, Liedtke C et al. (2012). Melanoma antigen family A identified by the bimodality index defines a subset of triple negative breast cancers as candidates for immune response augmentation. Eur J Cancer;48(1):12-23. PMID: 21741824
	\item Popovici V, Chen W, Gallas BG, Hatzis C et al. (2010). Effect of training-sample size and classification difficulty on the accuracy of genomic predictors. Breast Cancer Res;12(1):R5. PMID: 20064235
	\item Rody A, Karn T, Liedtke C, Pusztai L et al. (2011). A clinically relevant gene signature in triple negative and basal-like breast cancer. Breast Cancer Res, 6;13(5):R97. PMID: 21978456	
	\item Shi L, Campbell G, Jones WD, Campagne F et al. (2010). The MicroArray Quality Control (MAQC)-II study of common practices for the development and validation of microarray-based predictive models. Nat Biotechnol;28(8):827-38. PMID: 20676074
}
}

\keyword{datasets}

