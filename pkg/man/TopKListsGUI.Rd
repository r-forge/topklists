\name{TopKListsGUI}
\alias{TopKListsGUI}
\title{
TopKListsGUI for inference and visualization
}
\description{
This function opens a RGUI window and allows the user to select the parameters for the distance \eqn{\delta}, for the pilot sample size \eqn{\nu}, and the threshold for the \code{aggmap} presentation. Based on the selected parameters, a pairwise estimation of the top-\emph{k} lists overlap is performed with switching reference lists (e.g. L1 with L2 and L2 with L1). Based on these estimates, all involved lists are truncated and the maximum of the \eqn{\hat{k}_i} estimates is selected by default. The individual results of each combination of ranked lists are displayed automatically in the RGUI window and saved to a prespecified folder. After the maximal truncation point is estimated, \code{CEMC} algorithm is applied to generate the final list of top-\emph{k} objects.

The consolidated top-\emph{k} list results can be displayed in three different formats controlled by tabs: 

	(1)'Aggregation map': It is a special type of heatmap, where the truncated lists are ordered from left to right, from the one with the largest overlap with all the others to the one with the lowest overlap. For details see the description of the \code{aggmap}.
	(2) 'Summary table': An interactive table that displays the set of overlapping objects in dependence of the selected parameters. 
	(3) 'Venn-diagram & Venn-table': A Venn-diagram and a Venn-table of the overlapping objects based on all truncated input lists. 
}
\usage{
TopKListsGUI(lists, autorange.delta = FALSE, override.errors = TRUE,
venndiag.pdf.size = c(7, 7), venndiag.size = c(380, 420),
gui.size = c(900, 810), directory = NULL, venndiag.res = 70, aggmap.res = 100)
}
\arguments{
  \item{lists}{A data frame of ordered lists of objects - rows represent the objects, columns represent the individual input lists
}
\item{autorange.delta}{If TRUE, results for all \eqn{\delta} values leading to the  \eqn{\hat{k}_i} estimates (as function of the parameters specified by the user) will be displayed (note, for a large list length \emph{N} this can take quite some time, because the \eqn{\delta} range examined is from 0 to \code{nrow(lists)} by steps of 1). If FALSE, then the user defines the \eqn{\delta} range (min and max) in the RGUI window 
}
\item{override.errors}{If TRUE, errors due to an inappropriate value selection are overridden and the calculation continues for all other \eqn{\delta} values, but there is no output for viewing in the GUI corresponding to that combination of values which has caused the error; defaults to TRUE 
}
\item{venndiag.pdf.size}{A numeric vector defining the width and height of the Venn-diagram plot - it is passed to the \code{pdf()} function that saves the plot; defaults to c(7, 7)
}
\item{venndiag.size}{A numeric vector defining the width and height of the Venn-diagram plot - it is passed to the \code{png()} function that saves the plot; defaults to c(380, 420)
}
\item{gui.size}{A numeric vector defining the width and height of the RGUI to be displayed; defaults to c(900, 810)
}
\item{directory}{Specification of the name of the directory where the results and plots should be saved (including some temporary files required for the calculations). If kept \code{NULL} a directory called "TopKLists-temp" will be created in \code{getwd()}. In that case a write permission is needed for \code{getwd()}.
}
\item{venndiag.res}{A number defining the resolution of the Venn-diagram plot - this argument is passed to the \code{res = argument} of the \code{png()} function saving the diagram, but affects also the RGUI display; defaults to 70
}
\item{aggmap.res}{A number defining the resolution of the \code{aggmap} plot - this argument is passed to the \code{res = argument} of the \code{png()} function saving the plot, but affects also the RGUI display; defaults to 100
}
}
\value{
RGUI window with three tabs:

'Aggregation map':  For an index \eqn{p=1,2,\ldots$, $L-1} aggregation levels (groupings of top lists) are combined in one display. For each group of \emph{L-p} truncated lists down to the smallest group consisting of just one pair of lists, (1) an arbitrary reference list ("ground truth") is selected under the condition that it comprises \eqn{\max_i(\hat{k}_i)} items among all pairwise comparisons in the group of rankings, (2) symbols of its \eqn{\max_i(\hat{k}_i)} items are printed vertically from the highest to the lowest rank position, and (3) the aggregation information for all remaining \emph{L-p} rankings in the group is added, ordered according to descending list length. 

'Summary table': An interactive table that displays all overlapping (grey) objects based on the truncated list comparison. Rank sum per object and frequency of each object in the input lists or truncated lists are calculated over all compared lists. The first column denotes if an object was selected by the CEMC algorithm for the final set of common objects. The table can be ordered according to any of the displayed columns. 

'Venn-diagram & Venn-table': The Venn-diagram and the Venn-table display the rank intersection of the identified top-\emph{k} objects in two different formats. 

These tabs automatically save all plots and tables into the specified directory.

The following additional exploratory features are implemented:

'Deltaplot' (see \code{deltaplot}): For a preselected range of \eqn{\delta}'s and all list pairs, an exploratory plot of rank discordance is created and saved (function not part of the RGUI window).

'Mdelta' (see \code{deltaplot}): For a preselected range of \eqn{\delta}'s and all list pairs, Delta-matrices are created and saved (function not part of the RGUI window) in one \code{rdata} object (\code{Mdelta.rdata}). Each delta-matrix is saved individually in a tab delimited \code{.txt}-file.
}
\author{
Eva Budinska <budinska@iba.muni.cz>, Karl G. Kugler <kg.kugler@gmail.com>, Michael G. Schimek <michael.schimek@medunigraz.at>
}
\seealso{
  \code{\link{CEMC}}
}
\examples{
##using a GUI
\dontrun{
data(TopKGUISampleInput)
TopKListsGUI(lists)
# Run for nu=10, delta from 2 to 10
}

}

