\name{TopKListsGUI}
\alias{TopKListsGUI}
\title{
Visualizations of TopKListsGUI Aggmap
}
\description{
This function opens an RGUI window and allows the user to set the variable values (nu, delta range, lists, threshold) and perform a pairwise estimation of the top K list overlap in both directions (e.g. L1 with L2 and L2 with L1). Based on this estimate, lists are truncated at a the maximum k_hat estimate and compared. 
The result for each combination of variables is displayed automatically in the RGUI window and saved in the selected folder. 
The top k list overlap results are presented in three tabs: 
'Aggregation map' tab contains a special heatmap, where truncated lists are ordered from left to right, from the one with the largest overlap with others to the one with the lowest overlap (e.g. in an example of three lists if L2 vs L1 gives the largest estimate of top k, the lists will be ordered L2, L1, L3 and than if L1 vs L3 gives the largest k estimate than L3 vs L1 comparison, the following list order will be L1 vs L3).  
'Summary table': An intereactive table that shows all overlapping objects in the truncated list comparison. 
'Venn-diagram & Venn-table': This tab shows a Venn diagram and lists overlapping objects between lists. 
}
\usage{
TopKListsGUI(lists, autorange.delta = FALSE, override.errors = FALSE, aggmap.size = c(870, 440), aggmap.res = 100, venndiag.size = c(380, 420), venndiag.res = 70, gui.size = c(900, 810), directory = "TopKListsGUI_temp")
}
\arguments{
  \item{lists}{A data frame of ordered lists of genes - rows represent genes, columns represent individual lists
}
\item{autorange.delta}{If TRUE, results for all delta values leading to an k_max estimate (as a function of other variables specified by the user) will be displayed (attention, if gene lists are long, this can take long time, as the delta range examined is from 0 to nrow(lists) by 1) if FALSE, user defines the delta range (min and max) in the RGUI window. 
}
\item{override.errors}{If TRUE, errors due to inappropriate value selection are overridden and the calculation continues for all delta values, but there is no data to view in the GUI displayed for the combination of values causing the error. Defaults to FALSE. 
}
\item{aggmap.size}{A numeric vector defining the width and height of the aggmap plot - it is passed to the png() function that saves the plot. Defaults to c(870, 440).
}
\item{aggmap.res}{A number defining the resolution of the aggmap plot - This argument is passed to the res = argument of the png() function saving the plot, but affects also the RGUI display. Defaults to 100.
}
\item{venndiag.size}{A numeric vector defining the width and height of the venn diagram plot - it is passed to the png() function that saves the plot. Defaults to c(380, 420).
}
\item{venndiag.res}{A number defining the resolution of the venn diagram plot - This argument is passed to the res = argument of the png() function saving the diagram, but affects also the RGUI display. Defaults to 70.
}
\item{gui.size}{A numeric vector defining the width and height of the RGUI to be displayed. Defaults to c(900, 810).
}
\item{directory}{Specify the name of the directory where the results and plots should be saved. If the directory does not exist, it creates one. Default name is "TopKListsGUI_temp".
}
}
\value{
RGUI window with three tabs:

'Aggregation map':  If a gene in the truncated leading list on the left (e.g. L2) was present in the truncated list it is compared to (e.g. L1 or L3), it is colored grey. If not, it is colored white. Lower triangle is colored in function of the rank distance of the gene in the two compared lists. 

'Summary table': An intereactive table that shows all overlapping (grey) objects in the truncated list comparison. Rank sum per object and frequency of each object in the input lists or truncated lists is calculated over all the lists compared. Table can be ordered according to any of the columns displayed. 

'Venn-diagram & Venn-table': This tab shows a Venn diagram of overlapping objects from summary table and a table of overlaps between all lists combinations. 

This function also automatically saves all plots and tables into the selected directory.

}
\author{
Eva Budinska
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
