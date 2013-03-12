Borda <- function(input,space){
#Compute Borda scores and ranks based on four aggregation function
        nList=length(input)
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

        rank=matrix(0,nrow=N,ncol=nList)

        for (i in 1:N){
                #Find the name of the ith element in L
                a=L[i]
                #Find the rank in each list or set it to NA
                for (l in 1:nList){
                        present=sum(space[[l]]==a)
                        if (present==1){ #element found
                                index=rev(order(space[[l]]==a))[1]
                                rank[i,l]=e.Topk[[l]][index,2]
                                }
                        if (present==0) #element not found
                                rank[i,l]=NA
                        }
                }

        #Aggregate results
        aggreg.function=c("mean","median","geomean","l2norm")
        allranks=data.frame(matrix(0,nrow=N,ncol=length(aggreg.function)))
        names(allranks)=aggreg.function
        allfun=data.frame(matrix(0,nrow=N,ncol=length(aggreg.function)))
        names(allfun)=aggreg.function
        for (k in 1:length(aggreg.function)){
            allfun[,k]=sort(apply(rank,1,aggreg.function[k], na.rm = TRUE))
            allranks[,k]=L[order(apply(rank,1,aggreg.function[k], na.rm = TRUE))]
                }
        results=list(allranks,allfun)
        names(results)=c("RankedObjects", "BordaScores")
        return(results)
        }

geomean <- function(x, na.rm=TRUE){
        return(exp(mean(log(x),na.rm=TRUE)))
        }

KendallCriterion <- 
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
                        #not in one of the spaces
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

l2norm <- function(x, na.rm=TRUE){
        return(mean(abs(x)^2,na.rm=TRUE))
        }

MC.ranks <- function(elements, trans, a=0.15, delta=10^-15){
#Compute rankings based on the transition matrix from a MC algorithm
        n=nrow(trans)
        trans=trans*(1-a)+a/n
        A=matrix(0,nrow=n,ncol=n)
        for (i in 1:n) A[i,i]=1
        diff=1
        count=0
        while (diff >delta){
                A1=A%*%trans
                diff=max(A1-A)
                A=A1
                count=count+1
                }
        results=list(count, rev(sort(A[1,])), elements[rev(order(A[1,]))])
        names(results)=c("Iterations", "StationaryProbabilities", "RankedObjects")
        return(results)
        }

trans.matrix <- function(input, space){
#Compute transition matrices for all three MC algorithms
        #Both input and space are lists of the same length=nList
        nList=length(input)
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
        MC1=MC2=MC3=matrix(0,nrow=N, ncol=N)

        #Build lookup table
        lookup=matrix(c(rep(L,rep(N,N)),rep(L,N)),ncol=2)
        lookup=lookup[lookup[,1]!=lookup[,2],]

        #Build MC matrix
        for (i in 1:(nrow(lookup))){
                a=lookup[i,1]
                b=lookup[i,2]
                found=0; nn=0;
                for (l in 1:nList){
                        present=sum(space[[l]]==a)+sum(space[[l]]==b)
                        if (present==2){ #both elements in list
                        found=found+1
                        index1=rev(order(space[[l]]==a))[1]
                        index2=rev(order(space[[l]]==b))[1]
                        nn=nn+sum(e.Topk[[l]][index1,2]>e.Topk[[l]][index2,2])
                        }
                        }
                index1=rev(order(L==a))[1]
                index2=rev(order(L==b))[1]
                if (found!=0){MC1[index1,index2]=ceiling(nn/found)
                              MC2[index1,index2]=floor(nn/found+0.5)
                              MC3[index1,index2]=nn/found
                             }
                if (found==0) MC1[index1,index2]=MC2[index1,index2]=MC3[index1,index2]=0.5
                }
        MC1=MC1/N
        MC2=MC2/N
        MC3=MC3/N
        for (i in 1:N){
                MC1[i,i]=1-sum(MC1[i,-i])
                MC2[i,i]=1-sum(MC2[i,-i])
                MC3[i,i]=1-sum(MC3[i,-i])
                }
        return(list(e.Topk,L,N,MC1,MC2,MC3))
}

