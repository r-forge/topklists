\name{kendallModified}
\alias{kendallModified}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
kendallModified
}
\description{
TBD
}
\usage{
kendallModified(input, aggregate, space = NULL, p = 0.5, w = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{
%%     ~~Describe \code{input} here~~
}
  \item{aggregate}{
%%     ~~Describe \code{aggregate} here~~
}
  \item{space}{
%%     ~~Describe \code{space} here~~
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
function (input, aggregate, space = NULL, p = 0.5, w = NULL) 
{
    if (missing(input)) 
        stop("You need to input the individual top-k lists")
    if (missing(aggregate)) 
        stop("You need to have the aggregate top-k list")
    if (is.null(space) == TRUE) {
        common = sort(unique(unlist(input)))
        space = vector("list", length(input))
        for (i in 1:length(input)) space[[i]] = common
    }
    nList = length(input)
    if (is.null(w) == TRUE) 
        w = rep(1, nList)
    e.Topk = space
    for (l in 1:nList) {
        n = length(input[[l]])
        e.Topkl = rep(n + 1, length(space[[l]]))
        e.Topkl[match(input[[l]], space[[l]])] = 1:n
        e.Topk[[l]] = e.Topkl
    }
    L = unique(unlist(input))
    N = length(L)
    e.aggregate = matrix(0, nrow = N, ncol = 2)
    e.aggregate[, 1] = L
    e.aggregate[, 2] = length(aggregate) + 1
    e.aggregate[match(aggregate, e.aggregate), 2] = 1:length(aggregate)
    dall = 0
    for (l in 1:nList) {
        newSpace = unique(c(space[[l]], e.aggregate[, 1]))
        M = length(newSpace)
        scale = M * (M - 1)/2
        e.newSpace = as.data.frame(matrix(0, nrow = M, ncol = 3))
        e.newSpace[, 1] = newSpace
        e.newSpace[, 2:3] = 99999
        e.newSpace[match(space[[l]], newSpace), 2] = e.Topk[[l]]
        e.newSpace[match(e.aggregate[, 1], newSpace), 3] = e.aggregate[, 
            2]
        n = nrow(e.newSpace)
        long = rep(as.numeric(e.newSpace[1:(n - 1), 2]), (n - 
            1):1)
        allranks = matrix(0, nrow = length(long), ncol = 4)
        allranks[, 1] = long
        allranks[, 3] = rep(as.numeric(e.newSpace[1:(n - 1), 
            3]), (n - 1):1)
        for (i in 1:(n - 1)) {
            allranks[((i - 1) * n - i * (i - 1)/2 + 1):(i * n - 
                i * (i + 1)/2), 2] = as.numeric(e.newSpace[(i + 
                1):n, 2])
            allranks[((i - 1) * n - i * (i - 1)/2 + 1):(i * n - 
                i * (i + 1)/2), 4] = as.numeric(e.newSpace[(i + 
                1):n, 3])
        }
        notmissing = (allranks[, 1] != 99999) * (allranks[, 2] != 
            99999) * (allranks[, 3] != 99999) * (allranks[, 4] != 
            99999)
        dall = dall + (sum(((allranks[, 1] - allranks[, 2]) * 
            (allranks[, 3] - (allranks[, 4])) < 0) * notmissing + 
            ((allranks[, 1] - allranks[, 2]) * (allranks[, 3] - 
                (allranks[, 4])) == 0) * notmissing * p + (1 - 
            notmissing) * p))/scale
    }
    names(dall) = "Modified Kendall Distance"
    return(dall)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
