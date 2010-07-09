\name{ExcCEMC}
\alias{ExcCEMC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Execute CEMC
}
\description{
Execute CEMC by providing information in a text file
}
\usage{
ExcCEMC(para = "para.txt")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{para}{
text file containing parameters
}
}
\details{
There are two types of input files.Note that in the input files, anything after # denote comments. For information on the format, check the Vignette 'ExcCEMC'.

}
\value{
Results are either written to console or in a file (See Details section)
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Shili Lin <shili@stat.ohio-state.edu>
}
\note{
For samples of input txt files, see inst directory of this package
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## to demonstrate the use we have to temporarilly change the working directory
mypath <- system.file(package="TopKLists")
temp <- getwd()

setwd(paste(mypath, "extdata", sep="/"))

# for parameters of call see "para.txt" in the library "extdata" dir
ExcCEMC(para="para.txt")

##set wd back to normal
setwd(temp)
}
\keyword{ robust }
