<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.0 Transitional//EN">
<HTML>
<HEAD>
	<META HTTP-EQUIV="CONTENT-TYPE" CONTENT="text/html; charset=utf-8">
	<TITLE>TopKLists</TITLE>
	<link rel="stylesheet" type="text/css" href="incl/main.css">
</HEAD>
<BODY LANG="en-US" DIR="LTR">
<div style="width: 800px; text-align: center">
<H1>Welcome to the TopKLists page!</H1>
<br><br> 
<img src="img/website_pic.png" alt="Welcome to the TopKLists website"  width=800px />
</div>
<H2>Details</H2>
<div style="max-width:750px">
<P>For multiple ranked input lists (full or partial) representing the same set
of N objects, the package <strong>Top<font color="FireBrick">K</font>Lists</strong> offers 
</P>
<OL>
	<LI><P>statistical inference on the lengths of informative top-k
	lists</P>
	<LI><P>stochastic aggregation of full or partial lists</P>
	<LI><P>graphical tools for the statistical exploration of input
	lists, and for the visualization of aggregation results</P>
</OL>
Results that can be obtained are a sublist of highly conforming objects in a new rank order, an aggregation map of the top-k objects, or a Venn-diagram of the top-k objects. Typical applications are rank-based data integration and meta-analyses.
<P>The <B>latest stable build</B> can be found on <A HREF="http://cran.r-project.org/web/packages/TopKLists/index.html" TARGET="_blank">CRAN</A>.</P>
<P>You can find the <STRONG>latest news</STRONG> and the <STRONG>project summary page</STRONG> <STRONG><A HREF="https://r-forge.r-project.org/projects/topklists/" TARGET="_blank"><SPAN STYLE="text-decoration: none"><SPAN STYLE="font-weight: normal">here</SPAN></SPAN></A></STRONG>.
</P>
<P>For a <STRONG>detailed description</STRONG> of the functions provided by the package see the
<A HREF="http://cran.r-project.org/web/packages/TopKLists/TopKLists.pdf" TARGET="_blank">reference manual</A>.</P>
<P>When using TopKLists for the first time it is helpful to make use of <B>reproducible examples</B>, which can be found in the <A HREF="http://cran.r-project.org/web/packages/TopKLists/vignettes/TopKLists.pdf" TARGET="_blank">Vignette</A> 
(comprises the modules TopKInference, TopKSpace, TopKGraphics;  the graphical user interface TopKListsGUI is described as well).</P>
</div>
<H2>Show case</H2>
<div style="max-width:750px">
<table style="border: 1">
<tr><td>Comparison of miRNA measurements</td><td><a href="showcase_miRNA/topklists-miRNA.html">website</a></td><td><a href="showcase_miRNA/topklists-miRNA.R">R code</a></td><td><a href="showcase_miRNA/topklists-miRNA-sampledata.zip">data</a></td></tr>
</table>
</div>
<H2>Citation Information</H2>
<div style="max-width:750px">
<P>When using TopKLists in your work please cite:
<UL>
	<li><a href="http://www.degruyter.com/doi/10.1515/sagmb-2014-0093" target="_blank">Schimek, M. G. et al. (2015). TopKLists: a comprehensive R package for statistical inference, stochastic aggregation, and visualization of multiple omics ranked lists. Stat Appl Genet Mol Biol. 2015 Jun;14(3):311-6. doi: 10.1515/sagmb-2014-0093.</a></li>
</UL>
When referring to the module TopKInference in your work please cite:
<UL>
	<li>Hall, P. and Schimek, M. G. (2012).  Moderate deviation-based inference for random degeneration in paired rank lists. J. Amer. Statist. Assoc., 107, 661-672.</li>
</UL>
When referring to the module TopKSpace in your work please cite:
<UL>
	<li>Lin, S. and Ding, J. (2009). Integration of ranked lists via Cross Entropy Monte Carlo with applications to mRNA and microRNA studies. Biometrics, 65, 9-18.</li>
</UL>
</div>
<H2>Examples for the Usage</H2>
<div style="max-width:750px">
<UL>
  <li>Schimek, M. G., Budinska, E., Kugler, K. and Lin, S. (2011). Package "TopKLists" for rank-based genomic data integration.
    Proceedings of CompBio 2011, 434-440, DOI: 10.2316/P.2011.742-032.</li>
  <li>Schimek, M. G. and Bloice, M. (2012). Modelling the rank order of Web search engine results.
	In Komarek, A. and Nagy, S. (eds). Proceedings of the 27th International Workshop on Statistical Modelling. (e-book ISBN 978-80-263-0250-6), Vol. 1, 303-308. </li>
  <li>Schimek, M. G., Mysickova, A. and Budinska, E. (2012). An inference and integration approach for the consolidation of ranked lists.
	Communications in Statistics - Simulation and Computation, 41:7, 1152-1166. 
	</li>
</UL>
</div>
<H2>Other Links</H2>
<UL>
	<LI><P><A HREF="http://www.stat.osu.edu/~statgen/SOFTWARE/TopKCEMC/" TARGET="_blank">TopKCEMC
	and TopKSpace webpages of Shili Lin</A></P>
</UL>
</BODY>
</HTML>
