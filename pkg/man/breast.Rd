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
 \item Pusztai: \href{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE20271}{GEP GSE20271}
}
}
\references{

	TransBig: Desmedt C, Piette F, Loi S, Wang Y et al. (2007).  Strong time dependence of the 76-gene prognostic signature for node-negative breast cancer patients in the TRANSBIG multicenter independent validation series. Clin Cancer Res, 1;13(11):3207-14. PMID: 17545524
	
	Pusztai: Tabchy A., Valero V., Vidaurre T., Lluch A. et al. (2010). Evaluation of a 30-gene paclitaxel, fluorouracil, doxorubicin, and cyclophosphamide chemotherapy response predictor in a multicenter randomized trial in breast cancer. Clin Cancer Res 1;16(21):5351-61. PMID: 20829329
	
	MDACC: Shi L, Campbell G, Jones WD, Campagne F et al. (2010). The MicroArray Quality Control (MAQC)-II study of common practices for the development and validation of microarray-based predictive models. Nat Biotechnol;28(8):827-38. PMID: 20676074

}
\keyword{datasets}

