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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(input,space,aggregate,p,w){
#Compute Kenall's criterion
        nList=length(input)
        scale=rep(0,nList)
        for (l in 1:nList)
                scale[l]=length(input[[l]])*(length(input[[l]])-1)/2
        e.Topk=space #create extended topk lists of same length as corresp space
        for (l in 1:nList){
                temp=matrix(c(space[[l]],rep(0,length(space[[l]]))),ncol=2)
                n=length(input[[l]])
                for (i in 1:n){
                        index=rev(order(input[[l]][i]==temp[,1]))[1]
                        temp[index,2]=i
                        }
                temp[temp[,2]==0,2]=n+1
                e.Topk[[l]]=temp
                }
        vec=input[[1]]
        for (l in 2:nList) vec=c(vec,input[[l]])
        L=unique(vec)
        N=length(L)
        #create extended aggregate list of same length as union of topk lists
        e.aggregate=matrix(0,nrow=N,ncol=2)
        e.aggregate[,1]=L
        for (i in 1:(length(aggregate))){
                index=rev(order(aggregate[i]==e.aggregate[,1]))[1]
                e.aggregate[index,2]=i
                }
        e.aggregate[e.aggregate[,2]==0,2]=length(aggregate)+1
        dall=0
        for (l in 1:nList){
                #find S_1 U S_2 where S_1 is e.Topk[[l]][,1] and S_2 is
                # e.aggregate[,1]
                newSpace=unique(c(e.Topk[[l]][,1],e.aggregate[,1]))
                M=length(newSpace)
                e.newSpace=matrix(0,nrow=M,ncol=3)
                e.newSpace[,1]=newSpace
                for (i in 1:(nrow(e.Topk[[l]]))){
                        index=rev(order(e.Topk[[l]][i,1]==e.newSpace[,1]))[1]
                        e.newSpace[index,2]=e.Topk[[l]][i,2]
                        }
                e.newSpace[e.newSpace[,2]==0,2]=99999
                for (i in 1:(nrow(e.aggregate))){
                        index=rev(order(e.aggregate[i,1]==e.newSpace[,1]))[1]
                        e.newSpace[index,3]=e.aggregate[i,2]
                        }
                e.newSpace[e.newSpace[,3]==0,3]=99999
                #Now we have a matrix with 3 columns
                #1st col is names of elements in the union of the two spaces
                #2nd col is the rank (including k+1 and NA) of input list
                #3rd col is the rank (including k+1 and NA) of aggregate list
                d=0
                for (u in 1:(M-1))
                for (v in (u+1):M){
                        done=0
                if (sum(c(e.newSpace[u,2]==99999,e.newSpace[u,3]==99999,e.newSpace[v,2]==99999,e.newSpace[v,3]==99999))>0){
                        #not in one of the space
                        d=d+p
                        done=1
                        }
                        if ((done==0) && (((e.newSpace[u,2]-e.newSpace[v,2])*(e.newSpace[u,3]-e.newSpace[v,3]))==0)){
                        #both elements not ranked in one list
                        d=d+p
                        done==1
                        }
                        if ((done==0) && (((e.newSpace[u,2]-e.newSpace[v,2])*(e.newSpace[u,3]-e.newSpace[v,3]))!=0))
                        #at least one is ranked in each list
                        d=d+(((e.newSpace[u,2]-e.newSpace[v,2])*(e.newSpace[u,3]-e.newSpace[v,3]))<0)
                        }
                d=w[l]*d/scale[l]
                dall=dall+d
                }
        return(dall)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
