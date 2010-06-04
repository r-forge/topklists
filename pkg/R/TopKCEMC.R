    ####################################
    #Input file structure
    #Only one column without column name
    #The first row is the # number of lists
    #The following is the #  of lists section
    #Each section begin with number of items, then followed by these items
    #Total # of row, Total items(in all lists) + the number of list+1
    #By Feiyou Qiu, 11.2007
    ##################################
    #For Combination inputFile, all lists in one file
    #InputFile=listInputFileNames[1]

.FormatComInput<-function(InputFile){
        #InputFile="G_ALL.txt"
        #InputFile="input_list.txt"
        list_tab=read.table(InputFile,header=F)[,1]
        listNum=as.numeric(as.character(list_tab[1]))
        ItemNum=rep(0,listNum)
        listVector=list()
        mat_tmp=list()
        col_name=c("Item")
        AllCol=c()
        for(i in 1:listNum){
           startpos=sum(ItemNum[1:i-1])+1+i
           ItemNum[i]=as.numeric(as.character(list_tab[startpos]))
           mat_tmp[i]=list(cbind(c(1:ItemNum[i]),as.vector(list_tab[(startpos+1):(startpos+ItemNum[i])])))
           listVector[i]=as.vector(list(list_tab[(startpos+1):(startpos+ItemNum[i])]))
           col_name=c(col_name,paste("Group",i,sep=''))
           AllCol=union(AllCol,mat_tmp[i][[1]][,2])
        }
        #TotalItem=sum(ItemNum)
        list_mat=matrix(c(rep(NA,length(AllCol)*(listNum+1))),ncol=listNum+1,nrow=length(AllCol))
        colnames(list_mat)=col_name
        list_mat[,1]=AllCol
        for(i in 1:listNum){
          for(j in 1:ItemNum[i]){
              item=mat_tmp[i][[1]][j,2]
              list_mat[list_mat[,1]==item,i+1]=mat_tmp[i][[1]][mat_tmp[i][[1]][,2]==item,1]
          }
        }
        return(list_mat)
    }

#For Indvidual input files(each list has one file)
#listInputFileNames=c("G1.txt","G2.txt","G3.txt")
.FormatIndiviInput<-function(listInputFileNames){
        listNum=length(listInputFileNames)
        #InputFile="G_ALL.txt"
        #InputFile="input_list.txt"

        ItemNum=rep(0,listNum)
        #listVector=list()
        mat_tmp=list()
        col_name=c("Item")
        AllCol=c()
        NamesV=c()
        for(i in 1:listNum){
           filename=listInputFileNames[i]
           list_tab=read.table(filename,header=F)[,1]
           NamesV=cbind(NamesV,unlist(strsplit(filename, "\\."))[1])
           #startpos=sum(ItemNum[1:i-1])+1+i
           ItemNum[i]=length(list_tab)
           mat_tmp[i]=list(cbind(c(1:ItemNum[i]),as.vector(list_tab)))
           #listVector[i]=as.vector(list(list_tab[(startpos+1):(startpos+ItemNum[i])]))
           col_name=c(col_name,paste("Group",i,sep=''))
           AllCol=union(AllCol,mat_tmp[i][[1]][,2])
        }
        #TotalItem=sum(ItemNum)
        list_mat=matrix(c(rep(NA,length(AllCol)*(listNum+1))),ncol=listNum+1,nrow=length(AllCol))
        colnames(list_mat)=col_name
        list_mat[,1]=AllCol
        for(i in 1:listNum){
          for(j in 1:ItemNum[i]){
              item=mat_tmp[i][[1]][j,2]
              list_mat[list_mat[,1]==item,i+1]=mat_tmp[i][[1]][mat_tmp[i][[1]][,2]==item,1]
          }
        }
        return(list(list_mat,NamesV))
    }

.ConverListMatTotopK<-function(list_mat,NamesV){
        Topk_list=list()
        #NamesV = c("Luo","Welsh", "Dhama", "True", "Singh")
        colnum=dim(list_mat)[2]-1
        #collength=c()
        #total=c()
        for(i in 1:colnum){
              col_t=list_mat[!is.na(list_mat[,i+1]),c(1,i+1)]
              row_idx=sort(as.numeric(col_t[,2]),decreasing=F,index.return=T)[2]$ix
              Topk_list[NamesV[i]]=list(col_t[row_idx,1])
              #total=union(total, col_t[,1])
              #collength=c(collength,length(Topk_list[[i]]))
        }
        #k=length(total)
       #Topk=structure(.Data=Topk_list,.Names = NamesV ,class = "data.frame",row.names=c(1:max(collength)))
        return(Topk_list)
    }

ExcCEMC <-function(para="para.txt"){

      # direct output to a file
      #sink("CEMC_output.txt", append=FALSE, split=TRUE)

      if(file.exists(para)){
        source(para)
      }else print("Check the parameter file!")
      if(!exists("InputStyle")|| !exists("FileNames"))
       print("Please set the parameter for InputStyle and FileNames")
      #Set the default value
      if(!exists("k"))  OutputType="2"
      if(!exists("k")) {
         print("Please the specify the K in parameter file")
         return;
      }
      if(!exists("dm")) dm="s"         #distance measure, "s" for spearman, "k" for kendall(p=0)
      if(!exists("N")) N=1000         #Number of samples each iterate
      if(!exists("N1")) N1=1           #Number of best samples retained each iterate
      if(!exists("rho")) rho=0.1        #proportion of samples used to estimated new probability matrix
      if(!exists("e1")) e1=0.1         #for p
      if(!exists("e2")) e2=2           #for y
      if(!exists("w")) w=0.5          #weight of the new probabilities estimate each iterate
      if(!exists("b")) b=0            #parameter used in .blur function
      if(!exists("init.m")) init.m="p"     #initialization method, see the function init.p for details
      if(!exists("init.w")) init.w=0       #probability matrix initialization:
      if(!exists("d.w")) d.w=NULL       #weights for distances from different input lists

      listInputFileNames=unlist(strsplit(FileNames, "\\,"))
      if(InputStyle==0){
         #For Combination inputFile, all lists in one file
         list_mat=.FormatComInput(listInputFileNames[1])
         NamesV= names(as.data.frame(list_mat))[-1]
       }
       if(InputStyle==1){
        #For Indvidual input files(each list has one file)
        res=.FormatIndiviInput(listInputFileNames)
        list_mat=res[[1]]
        NamesV=res[[2]]
      }
      if(k==0){
         knum=dim(list_mat)[1]
      }else   knum=k
      TopK=.ConverListMatTotopK(list_mat,NamesV)

#        re=list(topK=result.ori,p=p,
#       output.other=data.frame(diff.p=diff.p, uc=uc,dist.s=dist.s,
#       dist.k=dist.k, iter=iter, time=time.use),
#       input.list=topK.input,input.par=data.frame(k=k,dm=I(dm),N=N,N1=N1,
#       rho=rho,e1=e1,e2=e2,w=w,b=b,init.m=I(init.m),init.w=init.w))
      re=CEMC(TopK,knum,dm,N,N1,rho,e1,e2,w,b,init.m,init.w,d.w,input.par=NULL)
      ##out put
      foo=.output(re)
      if(OutputType==1||OutputType==2){
        #ps <- pipe(paste("enscript -o ",Outputfile,".ps"),"w")  #"enscript -o ",
        #ps <- pipe("enscript -o tempout.ps","w")
        #capture.output(.output(re), file=ps)
        #close(ps)
        cat(foo,file="out.txt", sep = "\n")
      }
      if(OutputType==0||OutputType==2)  cat(foo, sep = "\n")
      #.output(re)

   }

.output<-function(re){
      zz <- textConnection("foo", "w")
      sink(zz)
      #if(Out_topK==1){
      string=""
      for(i in 1:length(re$topK)){
            string=paste(i,re$topK[i],sep=":")
            print(string)
      }
      #}
      #print(string)
      if(prob_out==1) print(re$p)
      #if(Out_input.list==1) print(re$input.list)
      #if(Out_input.par==1) print(re$input.par)
      #if(Out_other==1) print(re$output.other)

      sink()
      close(zz)
      return(foo)
}


############################################
#By Ding Jie, 11.2007
##############################################
`CEMC` <-
function(topK,k,dm="s",N=NULL,N1=NULL,rho=0.1,
         e1=0.1,e2=1,w=0.5,b=0,init.m="p",init.w=0, d.w=NULL,input.par=NULL){

  #topK: a list of several top K lists, may have different length,
  #      one column for each list, ended with NA's
  #k: length of needed list
  #N: number of samples each iterate
  #N1: number of best samples retained each iterate
  #rho: proportion of samples used to estimated new probability matrix
  #e1: for p
  #e2: for y
  #w: weight of the new probabilities estimate each iterate
  #b: parameter used in .blur function
  #dm: distance measure, "s" for spearman, "k" for kendall(p=0)
  #init.p: initialization method, see the function init.p for details
  #init.w: probability matrix initialization:
  #      (1-init.w) * uniform +  init.w * estimated from input lists
  #d.w: weights for distances from different input lists
  #input.par: dataframe of input parameters

  if (!is.null(input.par)) {
    for (p.n in names(input.par))
      assign(p.n,input.par[1,p.n])
  }

  time.start <- proc.time()

  topK.input <- topK

  a <- length(topK) # number of input top K lists
  k.a <- numeric(a) # lengths of input lists
  for (i in 1:a)
    k.a[i] <- length(topK[[i]])
  item <- sort(unique(unlist(topK))) # all items in the input top K lists
  n <- length(item)
  for (i in 1:a)
    topK[[i]] <- match(topK[[i]],item) # recode items from 1 to n

  #set default values for N and N1
  if (is.null(N))
    N <- 10 * n
  if (is.null(N1))
    N1 <- round(0.1 * N)

  p <- init.p(topK,n,k,init.m,init.w)
  p2 <- p
  y <- NA
  samp2 <- TopKSample.c(.blur(p,b),N1)
  iter <- 0
  y.count <- 0 # counter used in stopping rule

  if (is.null(d.w))
    d.w <- rep(1,a) #equal weight as default

  repeat {
    iter <- iter + 1
    samp <- cbind(samp2,TopKSample.c(.blur(p,b),N-N1))
    dist <- numeric(N)
    for (i in 1:a) {
      if (dm=="s")
        dist <- dist + spearman(topK[[i]],samp,n) * d.w[i]
      else if (dm=="k")
        dist <- dist + kendall.c(topK[[i]],samp,n,0) * d.w[i]
      else stop("Invalid distance measure")
    }
    y2 <- sort(dist)[round(N*rho)]
    samp2 <- samp[,order(dist)[1:N1]]
    samp3 <- samp[,dist<=y2]
    n.samp3 <- dim(samp3)[2]
    for (i in 1:k)
      for (j in 1:n)
        p2[j,i] <- sum(samp3[i,]==j)/n.samp3
    p2[,k+1] <- 1 - apply(p2[,1:k],1,sum)
    if (is.na(y)) {
      y <- y2
      p <- p*(1-w)+p2*w
    }
    else if (abs(y-y2) < e2) y.count <- y.count + 1
    else y.count <- 0
    y <- y2
    p <- p*(1-w)+p2*w
    if (y.count >= 5) break
  }
  result <- samp[,order(dist)[1]]
  dimnames(p) <- list(item,1:(k+1))

  diff.p <- mean(abs(p[result,1:k] - diag(k)))
  #difference between sorted p and the identity matrix

  uc <- mean(1-p[cbind(result,1:k)])
  #mean uncertainty for the final list

  result.ori <- item[result] #result in original coding

  dist.s <- 0
  dist.k <- 0
  for (i in 1:a) {
    dist.s <- dist.s + spearman(topK[[i]],result,n)
    dist.k <- dist.k + kendall.c(topK[[i]],result,n,0)
  }

  time.end <- proc.time()
  time.use <- sum((time.end-time.start)[-3])

  #Modified by feiyou Qiu 12/11/2007, add return value
  re=list(topK=result.ori,p=p,
       output.other=data.frame(diff.p=diff.p, uc=uc,dist.s=dist.s,
       dist.k=dist.k, iter=iter, time=time.use),
       input.list=topK.input,input.par=data.frame(k=k,dm=I(dm),N=N,N1=N1,
       rho=rho,e1=e1,e2=e2,w=w,b=b,init.m=I(init.m),init.w=init.w))
  #print(re)
  return(re)
}



`TopKSample` <-
function(p,N) {
  #Sampler to generate N top K lists according to p
  #p: matrix of dimension n*(k+1), n is the number of items
  #N: the number of samples

  k <- dim(p)[2] - 1
  n <- dim(p)[1]
  item <- 1:n
  samp <- matrix(0,nrow=k,ncol=N)

  i <- 1
  fail <- F
  repeat {
    samp[1,i] <- which(rmultinom(1,1,p[,i])==1)
    for (j in 2:k) {
      remain <- item[-samp[1:(j-1),i]]
      if (sum(p[remain,j]) == 0) {
        fail <- T
        break
      }
      #samp[j,i] <- remain[multinomial(p[remain,j])]
      samp[j,i] <- remain[which(rmultinom(1,1,p[remain,j])==1)]
    }
    if (!fail) i <- i + 1
    if (i == N+1) break
  }
  samp
}

`TopKSample.c` <-
function(p,N) {

  #this version calls a C subroutine for faster sampling

  k <- dim(p)[2] - 1
  n <- dim(p)[1]
  samp <- matrix(0,nrow=k,ncol=N)
  seed <- round(runif(1,10000,3000000))

  samp[] <- .C("topksamplec",as.double(p),as.integer(k),as.integer(n),
               as.integer(N),samp=as.integer(samp),as.integer(seed))$samp

  samp

}

`.blur` <-
function(p,b) {

  #p: vector of probabilites
  #b: exchange proportions, see below

  n <- dim(p)[1]
  k <- dim(p)[2]
  p2 <- matrix(0,nrow=n,ncol=k)

  p2[,1] <- p[,1]*(1-b) + p[,2]*b

  for (i in 2:(k-1))
    p2[,i] <- p[,i-1]*b+p[,i]*(1-2*b)+p[,i+1]*b

  p2[,k] <- p[,k]*(1-b) + p[,k-1]*b
  p2
}

`comp.rank` <-
function(input.list) {

  #calculate composite rank
  #input.list: a list of topk lists, may have different length

  n.list <- length(input.list) #number of lists
  nn.list <- unlist(lapply(input.list, length)) #lengths of individual lists

  item <- unique(unlist(input.list))
  n.item <- length(item)

  score <- matrix(0,nrow=n.item, ncol=n.list)

  for (i in 1:n.list) {
    list.code <- match(input.list[[i]],item) #coded list
    score[,i] <- nn.list[i]+1 #scores for items not in the list
    score[list.code,i] <- 1:nn.list[i] #scores for item in the list
    score[,i] <- score[,i]/nn.list[i] #weighted score
  }

  c.score <- apply(score,1,sum) #composite scores
  item[order(c.score)] #order w.r.t the composite scores

}

`init.p` <-
function(topK,n,k,init.m="p",init.w=0) {

  #topK: a list of input lists, with items coded from 1 to n
  #n: total number of items
  #k: length of target list
  #init.m: initialization method, "p" point mass, "s" smooth
  #        "cp" point mass using composite ranks
  #        "cs" smooth using composite ranks
  #init.w: initialization weight

  p.u <- matrix(1/n,nrow=n,ncol=k+1)
  p.u[,k+1] <- 1-k/n

  a <- length(topK)

  if (init.m %in% c("p","s")) {
    p.e <- matrix(0,nrow=n,ncol=n)
    for (i in 1:a) {
      k.a <- length(topK[[i]])
      p.e[cbind(topK[[i]],1:k.a)] <- 1/a + p.e[cbind(topK[[i]],1:k.a)]
      p.e[(1:n)[-match(topK[[i]],1:n)],(k.a+1):n] <-
        1/((n-k.a)*a) + p.e[(1:n)[-match(topK[[i]],1:n)],(k.a+1):n]
    }
#    if ((k+1)<n)
#      p.e <- cbind(p.e[,1:k],apply(p.e[,(k+1):n],1,sum))
    }
  else if (init.m %in% c("cp","cs")) {
    p.e <- matrix(0,nrow=n,ncol=n)
    c.list <- comp.rank(topK)
    p.e[cbind(c.list,1:n)] <- 1
#    p.e[,(k+1):n] <- (1-apply(p.e[,1:k],1,sum))/(n-k)
    }
  else stop("Invalid initialization method")

  if (init.m %in% c("s","cs")) {
      dist <- 0
      for (i in 1:(a-1)) {
        for (j in (i+1):a) {
          dist <- dist + spearman(topK[[i]],topK[[j]],n)
        }
      }
      dist.avg <- dist/(a*(a-1)/2)/n #average rank difference
      normal.c <- sum(dnorm(1:n,mean=(n+1)/2,sd=dist.avg/2))
                                     #normalizing constant
      weight <- matrix(0,nrow=n,ncol=n)
      for (i in 1:n) {
        weight[,i] <- dnorm(1:n,mean=i,sd=dist.avg/2)/normal.c
        weight[i,i] <- weight[i,i]+(1-sum(weight[,i])) #make weights sum to 1
      }
      p.e <- p.e %*% weight
     }
  p.e <- cbind(p.e[,1:k],1-apply(p.e[,1:k],1,sum))
  p.u*(1-init.w) + p.e*init.w
}

`kendall` <-
function(a,b,n=0,p=0) {

  #Kendall's tau distance between top K lists
  #n: total number of objects, numbered from 1 to n
  #a: a single top k list
  #b: top k list(s) need to be compared to list a
  #p: distance added for tied pair (potential problem when p != 0)

  if (n==0) {  #the lists are not coded for 1 to n
    item <- unique(c(a, as.vector(b)))
    n <- length(item)
    a <- match(a, item)
    if (is.vector(b))
      b <- match(b,item)
    else
      b <- matrix(match(b,item),nrow=dim(b)[1])
  }

  k.a <- length(a)
  rank.a <- rep(k.a+1,n)
  rank.a[a] <- 1:k.a

  if (is.vector(b)) {
    k.b <- length(b)
    n.b <- 1
    rank.b <- rep(k.b+1,n)
    rank.b[b] <- 1:k.b
    }
  else {
    k.b <- dim(b)[1]
    n.b <- dim(b)[2] # number of lists in b
    rank.b <- matrix(k.b+1,nrow=n,ncol=n.b)
    rank.b[cbind(as.vector(b),rep(1:n.b,each=k.b))] <- 1:k.b
    }

  dist <- numeric(n.b)
  d.a <- rank.a - rep(rank.a,each=n)
  for (i in 1:n.b) {
    if (n.b == 1) d.b <- rank.b - rep(rank.b,each=n)
    else d.b <- rank.b[,i] - rep(rank.b[,i],each=n)
    count <- table(c(sign(d.a)*sign(d.b),-1,0,1)) #make sure all 3 values exist
    dist[i] <- (count["-1"]-1)*1/2 + (count["1"]-1)*0/2 + (count["0"]-n-1)*p/2
  }
  dist
}

`kendall.c` <-
function(a,b,n=0,p=0) {

  #version that calls a C subroutine to speed up calculation
  #Kendall's tau distance between top K lists
  #n: total number of objects, numbered from 1 to n
  #a: a single top k list
  #b: top k list(s) need to be compared to list a
  #p: distance added for tied pair (potential problem when p != 0)

  if (n==0) {  #the lists are not coded for 1 to n
    item <- unique(c(a, as.vector(b)))
    n <- length(item)
    a <- match(a, item)
    if (is.vector(b))
      b <- match(b,item)
    else
      b <- matrix(match(b,item),nrow=dim(b)[1])
  }

  k.a <- length(a)
  rank.a <- rep(k.a+1,n)
  rank.a[a] <- 1:k.a

  if (is.vector(b)) {
    k.b <- length(b)
    n.b <- 1
    rank.b <- rep(k.b+1,n)
    rank.b[b] <- 1:k.b
    }
  else {
    k.b <- dim(b)[1]
    n.b <- dim(b)[2] # number of lists in b
    rank.b <- matrix(k.b+1,nrow=n,ncol=n.b)
    rank.b[cbind(as.vector(b),rep(1:n.b,each=k.b))] <- 1:k.b
    }

  dist <- numeric(n.b)
  .C("kendallc",as.integer(rank.a),as.integer(rank.b),as.integer(n),
     as.integer(n.b),as.double(p),dist=as.double(dist))$dist
}

`spearman` <-
function(a,b,n=0) {

#Spearman's footrule distance between two top K lists
#n: total number of objects, numbered from 1 to n
#a: a single top k list
#b: top k list(s) need to be compared to list a
#   each column is for one list

  if (n==0) {  #the lists are not coded for 1 to n
    item <- unique(c(a, as.vector(b)))
    n <- length(item)
    a <- match(a, item)
    if (is.vector(b))
      b <- match(b,item)
    else
      b <- matrix(match(b,item),nrow=dim(b)[1])
  }

  k.a <- length(a)
  if (is.vector(b)) k.b <- length(b)
  else k.b <- dim(b)[1]
  k <- max(k.a,k.b)

  rank.a <- rep(k+1,n)
  rank.a[a] <- 1:k.a

  if (is.vector(b)) {
    n.b <- 1
    rank.b <- rep(k+1,n)
    rank.b[b] <- 1:k.b
    }
  else {
    n.b <- dim(b)[2] # number of lists in b
    rank.b <- matrix(k+1,nrow=n,ncol=n.b)
    rank.b[cbind(as.vector(b),rep(1:n.b,each=k.b))] <- 1:k.b
    }

  if (n.b == 1)
    d <- sum(abs(rank.a-rank.b)*((rank.a<=k.a) | (rank.b<=k.b)))
    # for each list in b, only calculate the distance using the items
    # exist in list a and/or the list in b
  else
    d <- apply(abs(rank.a-rank.b)*((rank.a<=k.a) | (rank.b<=k.b)),2,sum)
  d
}

