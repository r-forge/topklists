\name{KendallCriterion}
\alias{KendallCriterion}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute Kendall's criterion
}
\description{
Compute Kenall's criterion
}
\usage{
KendallCriterion(input, space, aggregate, p, w)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{
%%     ~~Describe \code{input} here~~
}
  \item{space}{
%%     ~~Describe \code{space} here~~
}
  \item{aggregate}{
%%     ~~Describe \code{aggregate} here~~
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
Jie Ding <ding@stat.osu.edu>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
L1=1:30
L2=c(1:10,31:40,11:15)
L3=c(1:10,16:20,11:15)
input=list(L1,L2,L3)
space=list(1:40,1:40,1:40)
bb1=Borda(input,space)
kc.ARM=KendallCriterion(input, space, bb1[[1]][, 1], p = 0.5,
w = c(2/(30 * (30 - 1)), 2/(25 * (25 - 1)), 2/(20 * (20 -+ 1))))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
