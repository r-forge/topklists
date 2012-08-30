# method for j0 estimation
# the algorithm consists of an ordered sequence of "test stages" s1, s2....
# stage sk is associated with an integer Jsk, which when k is odd, is a potential lower bound to j0
# Idata - input data is a vector of 0's and 1's


compute.stream<-function(Idata, const=0.251, v, r=1.2)
{
if(sum(Idata, na.rm=T)==length(Idata)) 
{
return(list(j0_est=length(Idata), reason.break="Idata is identity", Js=NA, v=v))

}else
{
zv 		= .moderate.deviation(const, v)
Js 		= c()
k 		= 1
pj.plus 	= c()
pj.minus 	= c()
h		= 0
reason.break="NA"

###########
repeat
{
if (floor(k/2)<(k/2))
{
	if (k==1)
	{	h = 1;
		j=h-1;
		v.last=v[k]
		if (j< length(Idata))
		{
			repeat 
			{
			j = j+1
		
			if ((j+v[k]-1) >= length(Idata)) {#print("(j+v[k]-1) >= length(Idata)")
								    break}
	
			# computing pjplus
			pj.plus = .pjplus(Idata, v[k], j)
			#print(pj.plus)
			#print(Idata[j:(j+v[k]-1)])
			
			# testing
			if ((pj.plus-0.5)<=zv) {#print("pj.plus-0.5<=zv")
							break}
			}#end repeat

		Js[k] = j	
		Js	
		if(reason.break!="NA"){break;}				

		}# end if j<length(Idata)

	} else{
		h = Js[k-1]-trunc(r*v.last)+1;
		j = h-1;
		#print(paste("j=",j, sep=""))
		if (j< length(Idata))
		{
			repeat 
			{
			j = j+1
		
			if ((j+v.last-1) >= length(Idata)) {#print("(j+v.last-1) >= length(Idata)")
									break}
	
			# computing pjplus
			pj.plus = .pjplus(Idata, v.last, j)
			#print(pj.plus)
			#print(Idata[j:(j+v.last-1)])
			
			# testing
			if ((pj.plus-0.5)<=zv) {#print("pj.plus-0.5<=zv")
							break}
			}
		if(reason.break!="NA"){break;}
		Js[k] = j	
		Js					
		# breaking the repeat loop condition 2
		if ((k-3)>=1) {if (Js[k-2]==Js[k] & Js[k-1]==Js[k-3]) {reason.break="Js[k-2]==Js[k] & Js[k-1]==Js[k-3]";
							break}
			 		 }
		}# end if

	}#end else
} # end is.odd

if ((floor(k/2)==(k/2)))
{
	h = Js[k-1]+trunc(r*v.last)
	j = h

	if (j>1)
	{
		repeat
		{
		j = j-1

		# computing pjminus
		if ((j-v.last+1) <= 0){break}
		
		pj.minus 	= .pjminus(Idata, v.last, j)
		#print(pj.minus)
		#print(Idata[(j-v.last+1):j])
			
		# testing
		if ((pj.minus-0.5)>zv) {#print("pj.minus-0.5>zv")
						break}	
		}# end repeat

	Js[k] = j
	
	# breaking the repeat loop condition 3
	#if (Js[k]-trunc(r*v)+1<=1) {reason.break="Js[k]-r*v<=1"; break}
	 if (Js[k]-trunc(r*v.last)+1<=1) {  v.temp=v.last; 
							while((Js[k]-r*v.temp<=1) & v.temp>0){v.temp=v.temp-1}
							v.last=v.temp}
	 if (v.last==1){reason.break = "v converged to 1"; v[k+1]=v.last ;break; #print("selected v too small")
				}
	 if (v.last==0){reason.break = "v converged to 0"; v[k+1]=v.last ;break; #print("selected v too small")
				}


	# breaking the repeat loop condition 1
	if (is.infinite(Js[k-1]) & is.infinite(Js[k])) 
		{reason.break="is.infinite(Js[k-1]) & is.infinite(Js[k])"; break}

	} # end if
	if(reason.break!="NA"){break;}
}

k=k+1
v[k]=v.last
#print(paste("k=",k, sep=""))

}#end repeat
if(v.last<=1 | reason.break!="Js[k-2]==Js[k] & Js[k-1]==Js[k-3]"){j0_est=NA} else {if ((floor(k/2)<(k/2)) & k>2) {j0_est = ceiling(Js[k-2]+0.5*v[k-2]-1)}else if ((floor(k/2)==(k/2)) & k>1){j0_est = ceiling(Js[k-1]+0.5*v[k-1]-1)} else{j0_est=NA}}
return(list(j0_est=j0_est, reason.break=reason.break, Js=Js, v=v))
}# end if sum(Idata)==length(Idata)
}

#Method of M. G. Schimek extended to count K(j0_est) for each pair of several ranked lists. 
#Lists have same length and have no missing values.
#Inputs:
#Data - input matrix, each coloum presents one list
#delta - the maximal distance of the objects in ranking
#v - parameter for estimating j0
#outputs: maximal estimated J0 through all combinations of 2 lists

#Data = lists
#delta=10
#v=10

j0.multi<-function(Data,d,v) {
  maxK=0
  K = c()
  for (i in 1:ncol(Data)){
  for (j in 1:ncol(Data)){
    if (i!=j) {
	ID=prepareIdata(Data[,c(i,j)],d=d)
	J= compute.stream(ID$Idata,v=v)$j0_est
	K = rbind(K, cbind(names(Data)[i], names(Data)[j], v,J,d))
	}# end for if
      }# end for j
    }# end for i
K = data.frame(K, stringsAsFactors=F)
names(K) = c("list1", "list2", "v", "j0_est","delta")
if(sum(is.na(K$j0_est))<nrow(K)){maxK = max(as.numeric(K$j0_est), na.rm=T)}else{maxK=NA}
return(list(maxK=maxK,K=K))
}

# Zv moderate deviation

.moderate.deviation <- function(C, v)
{
zv = sqrt((C*log(v, base=10))/v)
return(zv)
}

.pjminus <- function(Idata, v, j)
{
pj.minus = 1/v*(sum(Idata[(j-v+1):j], na.rm = TRUE))
return(pj.minus)
}

.pjplus <- function(Idata, v, j)
{
pj.plus = (1/v)*(sum(Idata[j:(j+v-1)], na.rm = TRUE))
return(pj.plus)
}

# function to prepare Idata from several rankings of several assessors
# x 	- data matrix, where columns represents rank order of genes from two different assessors, rows represents genes (clones)
# delta - the maximal distance of the objects in ranking between assessors
# num.omit - the maximal number of ommited genes from the analysis 
#
# the result is a object of type "Idata", which is a list containing Idata and the information
# about the distance (delta).

prepareIdata <- function(x, d)
{
rank.diff = c(1:nrow(x))-match(x[,1],x[,2])
Idata = abs(rank.diff)<=d
return(list(Idata = Idata, d = d))
}

is.odd <- function(k){
  return(k %% 2 != 0)
}

is.even <- function(k){
  return(!is.odd(k))
}

## Heatmap solution to HS rank order graphics representation
## By Eva Budinska, 3.12.2010
## Input variables 
#K- number of lists
# N maximal number of genes
# delta=10 maximal distance between genes used for j0 estimation
# v - nu value used for j0 estimation
# res.temp - data frame, result of HS algorithm for all pairs of lists, 4 columns, first and second column contain lists names that were compared, third column selected nu and fourth column estimated j0
# lists data frame containing three orderesd lists of genes that were compared (the input to HS algorithm) - column names obligatory

aggmap<-function(lists, res.j0.temp, N, K, d, v) {
  require(grid)
  wid=1
  hei=1
  res.temp =  as.matrix(res.j0.temp$K)
  max.j0.est = res.j0.temp$maxK

  local.max.j0.est = rep(K, K)
  names(local.max.j0.est) = names(lists)
  temp2=c()

  for (i in c(1:K))
  {
  if(i>1){res.temp.temp = res.temp[-c(which(!is.na(match(res.temp[,1],names(temp2)))),which(!is.na(match(res.temp[,2],names(temp2))))),]}else{res.temp.temp=res.temp}
  temp = tapply(as.numeric(res.temp.temp[,4]), res.temp.temp[,1], function(x) max(x, na.rm=TRUE))
  which.max.temp = which.max(temp)
  temp2 = c(temp2, temp[which.max.temp])
  }

  local.max.j0.est[names(temp2)]=order(temp2, decreasing=T)
  first.order.local.max.j0.est = order(local.max.j0.est)

  red.white.green = colorRampPalette(c("red", "orange","yellow","white"))(N+1)
  grid.newpage()

  # defining basic layout with two rows and one column
  pushViewport(viewport(layout=grid.layout(2, 1, widths=c(1,1),heights=c(0.1,0.9))))

  ##first row - defining header (delta, nu and color key level) 
	  pushViewport(viewport(layout.pos.row=1, layout.pos.col=1))
  ### first column - delta and nu 
		  pushViewport(viewport(layout=grid.layout(1, 1, widths=c(1,1),heights=rep(1,2))))
		  grid.text(paste("delta = ",d,",  v =",v,sep=""), x=0.2, y=0.5)
		  upViewport()		

  ### second column - color key
		  pushViewport(viewport(layout=grid.layout(1, 2, widths=c(1,1),heights=rep(1,2)))) # defines two columns in the second column
			  pushViewport(viewport(layout.pos.row=1, layout.pos.col=2)) # fills the second column with the key
				  pushViewport(plotViewport(margins=c(1.5,1,1,1)))
				  grid.xaxis(at=seq(0, 1, length=11), gp=gpar(cex=0.5,tck=0.1, line=0.1))
				  grid.points(seq(0, 1, length=N), rep(1, N), pch=15, gp=gpar(col=red.white.green))
				  popViewport()
			  popViewport()
		  upViewport()
	  upViewport() #coming back to the basic level

  ## second row - plotting genes
	  pushViewport(viewport(layout.pos.row=2, layout.pos.col=1))

  ### defining the layout in second row
		  pushViewport(viewport(layout=grid.layout(max.j0.est+v+1, 1+factorial(K)/2+2, widths=rep(1,1+factorial(K)/2+2)*(max.j0.est+v+1),heights=rep(1,1+factorial(K)/2+2)*(max.j0.est+v+1))))

			  for (g in c(1:(max.j0.est+v)))
				  {
				  pushViewport(viewport(width=wid, height=hei, layout.pos.col=1, layout.pos.row=g+1))
				  grid.text(g, gp=gpar(fontsize=8, col="black"))
				  popViewport()
				  }# end for g
		  h.sizes=c(K:2)
		  l=1
		  lists.to.remove = ""

	  for (first.list.name in names(temp2))
		  {
		  # selecting the order of pairs of lists to the first one
		  res.temp.selection = res.temp[which(res.temp[,1]==first.list.name),]
		  res.temp.selection = res.temp.selection[which(is.na(match(res.temp.selection[,2],lists.to.remove))),]
			  if(is.matrix(res.temp.selection)){
			  ranked.list.names = c(first.list.name, res.temp.selection[order(res.temp.selection[,4], decreasing=TRUE),2])
			  lists.rank = rank(ranked.list.names)
			  order.local.max.j0.est = order(as.numeric(res.temp.selection[,4]), decreasing=TRUE)
			  }else{
			  ranked.list.names = c(first.list.name,res.temp.selection[2])
			  lists.rank = rank(ranked.list.names)
			  order.local.max.j0.est = 1
			  }

		  gene.names.to.plot = lists[,first.list.name][1:(temp2[first.list.name]+v)]


		  ### plotting gene names ###
		  pushViewport(viewport(width=wid, height=hei, layout.pos.col=l+1,layout.pos.row=1))
		  grid.text(ranked.list.names[1], gp=gpar(fontsize=8, col="black"))
		  popViewport()

		  for (g in c(1:(length(gene.names.to.plot))))
				  {
				  pushViewport(viewport(width=wid, height=hei, layout.pos.col=l+1, layout.pos.row=g+1))
				  grid.rect(gp=gpar(col="black"))
				  grid.text(gene.names.to.plot[g], gp=gpar(fontsize=8, col="black"))
				  popViewport()
				  }# end for g
		  
		  for (k in c(2:K))
		  {
		  pushViewport(viewport(width=wid, height=hei, layout.pos.col=l+k, layout.pos.row=1))
		  grid.text(ranked.list.names[k], gp=gpar(fontsize=8, col="black"))
		  popViewport()
		  n.genes.to.plot = as.numeric(res.temp[res.temp[,1]==first.list.name & res.temp[,2]==ranked.list.names[k],4])+v
			  for (g in c(1:n.genes.to.plot))
			  {
			  j = match(gene.names.to.plot[g],lists[,ranked.list.names[k]])
			  distance = abs(g-j)
				  if (!is.na(j) & j<(n.genes.to.plot+v)) # if gene is in the second list AND is in the top j0+v genes of the second list
					  {i="grey"}else{i="white"}
				  
				  pushViewport(viewport(width=wid, height=hei, layout.pos.col=l+k, layout.pos.row=g+1))
				  if (length(distance)>0)
					  {
					  grid.polygon(x=c(0,1,1),y=c(1,1,0), gp=gpar(col="black", fill=red.white.green[distance+1]))#upper right triangle
					  grid.text(distance,x=0.8, y=0.6, gp=gpar(cex=0.7))
					  }else
					  {
					  grid.polygon(x=c(0,1,1),y=c(1,1,0), gp=gpar(col="black", fill="white"))#upper right triangle
					  grid.text("NA",x=0.8, y=0.6, gp=gpar(cex=0.7))
					  }
				  grid.polygon(x=c(0,0,1),y=c(1,0,0), gp=gpar(col="black", fill=i))#bottom left triangle
				  popViewport()
			  }# end for g
		  }# end for k
	  
			  l=l+K
			  K=K-1
			  lists.to.remove = c(lists.to.remove, first.list.name)

		  }# end for first list name


			  upViewport()
		  
	  upViewport()
} # end of the function
