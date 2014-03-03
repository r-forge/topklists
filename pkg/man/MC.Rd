\name{MC}
\alias{MC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MC(input, space = NULL, K = NULL, a = 0.15, delta = 10^-15)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{
%%     ~~Describe \code{input} here~~
}
  \item{space}{
%%     ~~Describe \code{space} here~~
}
  \item{K}{
%%     ~~Describe \code{K} here~~
}
  \item{a}{
%%     ~~Describe \code{a} here~~
}
  \item{delta}{
%%     ~~Describe \code{delta} here~~
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
function (input, space = NULL, K = NULL, a = 0.15, delta = 10^-15) 
{
    if (missing(input)) 
        stop("You need to input the top-k lists to be aggregated")
    if (is.null(space) == TRUE) {
        common = sort(unique(unlist(input)))
        space = vector("list", length(input))
        for (i in 1:length(input)) space[[i]] = common
    }
    out.trans = trans.matrix(input, space)
    N = length(out.trans[[1]])
    MC1 = MC.ranks(out.trans[[1]], out.trans[[2]], a, delta)
    MC2 = MC.ranks(out.trans[[1]], out.trans[[3]], a, delta)
    MC3 = MC.ranks(out.trans[[1]], out.trans[[4]], a, delta)
    if (is.null(K) == TRUE) 
        K = N
    else {
        if (K > N) 
            K = N
    }
    results = list(MC1[[3]][1:K], MC1[[2]], MC2[[3]][1:K], MC2[[2]], 
        MC3[[3]][1:K], MC3[[2]])
    names(results) = c("MC1.TopK", "MC1.Prob", "MC2.TopK", "MC2.Prob", 
        "MC3.TopK", "MC3.Prob")
    return(results)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
