\name{TopKListsGUI}
\alias{TopKListsGUI}
\title{
TopKListsGUI for inference and visualization
}
\description{
This function opens a RGUI window and allows the user to select the parameters for the distance delta, for the pilot sample size nu, and the threshold in the aggmap presentation. Based on the selected parameters, a pairwise estimation of the top-K lists overlap is performed with switching reference lists (e.g. L1 with L2 and L2 with L1). Based on these estimates, all involved lists are truncated at a the maximum k_hat estimate and then combined. The result for each combination of variables is displayed automatically in the RGUI window and saved in the selected folder. 

The consolidated top-K list results can be displayed in three different formats: 
(1) 'Aggregation map': It is a special type of heatmap, where the truncated lists are ordered from left to right, from the one with the largest overlap with all the others to the one with the lowest overlap. For details see the description of the \code{aggmap}.
(2) 'Summary table': An interactive table that displays the set of overlapping objects depending on the selected parameters. 
(3) 'Venn-diagram & Venn-table': A Venn-diagram and a Venn-table of the overlapping objects based on all truncated input lists. 
}
\usage{
TopKListsGUI(lists, autorange.delta = FALSE, override.errors = FALSE, aggmap.pdf.size = c(9, 8), venndiag.pdf.size = c(7, 7), venndiag.size = c(380, 420), aggmap.size = c(870, 440), gui.size = c(900, 810), directory = tempdir(), venndiag.res = 70, aggmap.res = 100, maxd = 500)
}
\arguments{
  \item{lists}{A data frame of ordered lists of objects - rows represent the objects, columns represent the individual input lists
}
\item{autorange.delta}{If TRUE, results for all delta values leading to an k_max estimate (as a function of the parameters specified by the user) will be displayed (note, for a large list length N this can take quite long time, because the delta range examined is from 0 to nrow(lists) by 1) if FALSE, the user is defining the delta range (min and max) in the RGUI window. 
}
\item{override.errors}{If TRUE, errors due to inappropriate value selection are overridden and the calculation continues for all delta values, but there is no data to view in the GUI displayed for the combination of values causing the error. Defaults to FALSE. 
}
\item{aggmap.pdf.size}{A numeric vector defining the width and height of the aggmap plot - it is passed to the pdf() function that saves the plot. Defaults to c(9, 8).
}
\item{venndiag.pdf.size}{A numeric vector defining the width and height of the Venn-diagram plot - it is passed to the pdf() function that saves the plot. Defaults to c(7, 7).
}
\item{venndiag.size}{A numeric vector defining the width and height of the Venn-diagram plot - it is passed to the png() function that saves the plot. Defaults to c(380, 420).
}
\item{aggmap.size}{A numeric vector defining the width and height of the aggmap plot - it is passed to the png() function that saves the plot. Defaults to c(870, 440).
}
\item{gui.size}{A numeric vector defining the width and height of the RGUI to be displayed. Defaults to c(900, 810).
}
\item{directory}{Specification of the name of the directory where the results and plots should be saved (including some temporary files used for calculation). Uses the \code{tempdir()} by default.
}
\item{venndiag.res}{A number defining the resolution of the Venn-diagram plot - this argument is passed to the res = argument of the png() function saving the diagram, but affects also the RGUI display. Defaults to 70.
}
\item{aggmap.res}{A number defining the resolution of the aggmap plot - this argument is passed to the res = argument of the png() function saving the plot, but affects also the RGUI display. Defaults to 100.
}
\item{maxd}{for deltaplot - specification of the maximal delta-value that should be explored in the delta plot. Defaults to 500.
}
}
\value{
RGUI window with three tabs:

'Aggregation map':  For an index \eqn{p=1,2,\ldots$, $L-1} aggregation levels (groupings of top lists) are combined in one display. For each group of \emph{L-p} truncated lists down to the smallest group consisting of just one pair of lists, (1) an arbitrary reference list ("ground truth") is selected under the condition that it comprises \eqn{\max_i(\hat{k}_i)} items among all pairwise comparisons in the group of rankings, (2) symbols of its \eqn{\max_i(\hat{k}_i)} items are printed vertically from the highest to the lowest rank position, and (3) the aggregation information for all remaining \emph{L-p} rankings in the group is added, ordered according to descending list length. 

'Summary table': An intereactive table that displays all overlapping (gray) objects based on the truncated list comparison. Rank sum per object and frequency of each object in the input lists or truncated lists are calculated over all compared lists. The table can be ordered according to any of the displayed columns. 

'Venn-diagram & Venn-table': The Venn-diagram and the Venn-table display the rank intersection of the identified top-K objects in two different formats. 

This function also automatically saves all plots and tables into the selected directory.

'Delta-plot' - (see \code{deltaplot}) For a preselected maxd (delta ranges from zero to maxd in steps of one) and all list pairs, an exploratory plot of rank discordance is created and saved (function not part of the RGUI window)
'Mdelta' - (see \code{deltaplot}) For a preselected maxd (delta ranges from zero to maxd in steps of one) and all list pairs, Delta-matrices are created and saved (function not part of the RGUI window) in one rdata object (Mdelta.rdata) Each Delta-matrix is saved individually in a tab delimited .txt file.
}
\author{
Eva Budinska <budinska@iba.muni.cz>, Karl G. Kugler <karl@eigenlab.net>, Michael G. Schimek <michael.schimek@medunigraz.at>
}
\examples{
##using a GUI
\dontrun{
data(TopKGUISampleInput)
TopKListsGUI(lists)
# Run for nu=10, delta from 2 to 10
}

}
\keyword{ graphics }
