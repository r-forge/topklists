\name{plotKendall}
\alias{plotKendall}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotKendall(input, all.aggregates, space = NULL, algorithm = NULL, p = 0.5, w = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{
%%     ~~Describe \code{input} here~~
}
  \item{all.aggregates}{
%%     ~~Describe \code{all.aggregates} here~~
}
  \item{space}{
%%     ~~Describe \code{space} here~~
}
  \item{algorithm}{
%%     ~~Describe \code{algorithm} here~~
}
  \item{p}{
%%     ~~Describe \code{p} here~~
}
  \item{w}{
%%     ~~Describe \code{w} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (input, all.aggregates, space = NULL, algorithm = NULL, 
    p = 0.5, w = NULL) 
{
    if (missing(input)) 
        stop("You need to input the individual top-k lists")
    if (missing(all.aggregates)) 
        stop("You need to have the aggregate top-k lists for comparison")
    if (is.null(space) == TRUE) {
        common = sort(unique(unlist(input)))
        space = vector("list", length(input))
        for (i in 1:length(input)) space[[i]] = common
    }
    nList = length(input)
    if (is.null(w) == TRUE) 
        w = rep(1, nList)
    n = length(all.aggregates)
    if (is.null(algorithm) == TRUE) 
        algorithm = paste(rep("alg", n), 1:n)
    kd = rep(0, n)
    for (i in 1:n) kd[i] = Kendall(input, all.aggregates[[i]], 
        space, p, w)
    names(kd) = algorithm
    plot(1:n, kd, type = "o", xaxt = "n", xlab = "Algorithm", 
        ylab = "Modified Kendall Distance")
    axis(1, at = 1:n, labels = algorithm)
    kdl = list(kd)
    names(kdl) = "Modified Kendall Distance"
    return(kdl)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
