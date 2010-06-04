\name{MC.ranks}
\alias{MC.ranks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute rankings based on a MC algorithm
}
\description{
Compute rankings based on the transition matrix from a MC algorithm
}
\usage{
MC.ranks(elements, trans, a = 0.15, delta = 10^-15)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{elements}{
%%     ~~Describe \code{elements} here~~
}
  \item{trans}{
%%     ~~Describe \code{trans} here~~
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
Shili Lin <shili@stat.ohio-state.edu>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#get sample data
data(TopKSpaceSampleInput)

mm=trans.matrix(input,space) #platform dependent
MC1=MC.ranks(mm[[2]],mm[[4]])
MC2=MC.ranks(mm[[2]],mm[[5]])
MC3=MC.ranks(mm[[2]],mm[[6]])

mm1=trans.matrix(input,input) #top-k spaces
MC11=MC.ranks(mm1[[2]],mm1[[4]])
MC21=MC.ranks(mm1[[2]],mm1[[5]])
MC31=MC.ranks(mm1[[2]],mm1[[6]])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
