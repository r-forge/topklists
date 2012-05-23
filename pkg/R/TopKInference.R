                                                                     
                                                                     
                                                                     
                                             
# The first part of the module TopKInference provides exploratory nonparametric inference for the estimation of the top-k list length of paired rankings.

# Implemented in function compute.stream() is the moderate deviation-based inference procedure for random degeneration in paired rankings due to Hall & Schimek (2010). This exploratory (because it depends on the choice of the pilot sample size nu) non-parametric procedure gives an estimate of the point of degeneration j0, where k=j0-1 denotes the estimated length of the top list. It allows for various (usually unknown) types of rank irregularities and list lengths in the magnitude of thousands of objects. Input to this procedure is a 0/1-vector of indicator variables Ij (a data stream) representing concordance between the two rankings of interest. Ij is a Bernoulli random variable for which either independence or weak m-dependence is assumend (from simulations we know that dependence is negligible). The function prepareIdata() allows to calculate the data stream from two ranked lists representing the same N objects.

# Idata - data input is a vector of 0's and 1's.
# const - constant C of the moderate deviation bound for testing (C>0.25 required, default value is 0.251).
# v - pilot sample size nu (a smoothing parameter) depends on the data (integer values v>=1 required, see Hall & Schimek (2010) for hints on its choice).
# r - technical constant r connected to the "test stages" of the iterative algorithm (r>1 required, default value is 1.2).

# Output is the estimate j0 and otherwise the break condition (illegal input data or reason for non-convergence).

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

#------------------------------------------------------------------------------------------

# The second part of the module TopKInference provides a method that allows us to combine compute.stream() results when there are more than two rankings (multiple input lists) representing the same set of N objects (missing rank assignments cannot be handled in the current version but is an option in future ones). The method is described in Schimek et al. (2011). Function j0.multi() calculates an overall estimate of the index k based on all pairwise combinations of L input lists.

#Data - input matrix, each coloum presents one list.
#delta - the tolerated (predefined) distance delta between rank positions of the individual objects in the input lists. The parameter delta takes the integer values 0, 1, 2,.... The Delta-plot of the module TopKGraphics is helping the user to select delta in practice (for details see Schimek & Budinska, 2010). 
#v - pilot sample size nu (a smoothing parameter) depends on the data (integer values v>=1 required, see Hall & Schimek (2010) for hints on its choice).

# Output is a combined estimate j0 for the L input lists as a function of the individual j0 (their maximum is default). 

j0.multi<-function(Data,delta,v) {
  maxK=0
  K = c()
  for (i in 1:ncol(Data)){
  for (j in 1:ncol(Data)){
    if (i!=j) {
	ID=prepareIdata(Data[,c(i,j)],delta=delta)
	J= compute.stream(ID$Idata,v=v)$j0_est
	K = rbind(K, cbind(names(Data)[i], names(Data)[j], v,J,delta))
	}# end for if
      }# end for j
    }# end for i
K = data.frame(K, stringsAsFactors=F)
names(K) = c("list1", "list2", "v", "j0_est","delta")
maxK = max(as.numeric(K$j0_est), na.rm=T)
return(list(maxK=maxK,K=K))
}

#------------------------------------------------------------------------------------------

#Below are various functions required by compute.stream() and j0.multi()

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

#------------------------------------------------------------------------------------------

# Calculation of the Idata input from the ranked lists of two or more assessors.
# x - data matrix, where the columns represent the rankings of the objects and the rows the object names.
# delta - the tolerated (predefined) distance delta between rank positions of the individual objects in the input lists. The parameter delta takes the integer values 0, 1, 2,.... The Delta-plot of the module TopKGraphics can help the user to select the parameter delta in practice (for details see Schimek & Budinska, 2010). We suggest to calculate data matrices x for a range of plausible delta values.
# num.omit - the number phi of discarded objects in each test stage when such an adjustment is necessary (see Hall & Schimek, 2010). 

# The output is an object of type "Idata", which is a list containing Idata and the information
# about the distance delta.

######################## QUESTION FOR EVA #############################
#
########### DOES prepareIdata NOW WORK FOR BOTH SIMPLE INPUT OF 2 RANKINGS AND FOR MULTIPLE 
########### RANKINGS ???????????
########### PLEASE REMEMBER MY EMAIL CONCERNING PROBLEMS WITH A 2 LISTS EXAMPLE AND KARLS 
########### EXAMPLE IN THE DOCUMENTATION
#
##################################################################################

prepareIdata <- function(x, delta)
{
rank.diff = c(1:nrow(x))-match(x[,1],x[,2])
Idata = abs(rank.diff)<=delta
return(list(Idata = Idata, delta = delta))
}

is.odd <- function(k){
  return(k %% 2 != 0)
}

is.even <- function(k){
  return(!is.odd(k))
}


######################## COMMENT FOR EVA AND KARL #############################
#
########### delta SELECTION BEFORE RUNNING prepareIdata REQUIRES THE Delta-PLOT !!!!
########### EVA, PLEASE INCLUDE HERE A FUNCTION PRODUCING THE DELTA-PLOT
########### IN THE FUTURE THIS SHOULD BE PART OF TopKGraphics
#
##################################################################################

#------------------------------------------------------------------------------------------

# The function aggmap() belongs to the module TopKGraphics. For convenience of the user we provide it here because TopKGraphics is not released yet.

# The aggregation map as implemented in the function aggmap() is a heatmap-like graphical representation of the combined result of several full or truncated (top k) input lists. It has was introduced in Schimek & Budinska (2010) for visual inspection of the aggregation outcome of L input lists.

# N - maximal number of objects (for outcome to be viewed/plotted in one page, N should not exceed 100)
# K - number of input lists
######################## COMMENT FOR EVA AND KARL #############################
#
########### THE NUMBER OF LISTS SHOULD BE L THROUGHOUT AND NOT K BECAUSE K WOULD BE 
########### CONFUSED WITH "TOP K". PLEASE REPLACE ALSO IN THE CODE !!!!!!!!!!
#
##################################################################################
# delta - the tolerated (predefined) distance delta between rank positions of the individual objects in the input lists.
# v - pilot sample size nu (a smoothing parameter) 
# lists - data frame containing all ranked lists of the ordered objects subject to comparison (the column names are obligatory).
# res.j0.temp - data frame of the results from the inference procedure for all combinations of input lists. It consists of four columns: the first and second column contain the input list names that were compared, the third column the applied nu, and the fourth column the resulting j0 estimate.

# The output is a color plot.

aggmap<-function(lists, res.j0.temp, N, K, delta, v) {
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

  ##first row - defining the header (delta, nu and color key level) 
	  pushViewport(viewport(layout.pos.row=1, layout.pos.col=1))
  ### first column - delta and nu 
		  pushViewport(viewport(layout=grid.layout(1, 1, widths=c(1,1),heights=rep(1,2))))
		  grid.text(paste("delta = ",delta,",  v =",v,sep=""), x=0.2, y=0.5)
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

  ## second row - plotting of objects
	  pushViewport(viewport(layout.pos.row=2, layout.pos.col=1))

  ### defining the layout in the second row
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
		  # selecting the order of pairs of lists with respect to the first (reference) one
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


		  ### plotting of the object names ###
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
				  if (!is.na(j) & j<(n.genes.to.plot+v)) # if the object is in the second list AND is among the top j0+v objects of the second list
					  {i="grey"}else{i="white"}
				  
				  pushViewport(viewport(width=wid, height=hei, layout.pos.col=l+k, layout.pos.row=g+1))
				  if (length(distance)>0)
					  {
					  grid.polygon(x=c(0,1,1),y=c(1,1,0), gp=gpar(col="black", fill=red.white.green[distance+1])) #upper right triangle
					  grid.text(distance,x=0.8, y=0.6, gp=gpar(cex=0.7))
					  }else
					  {
					  grid.polygon(x=c(0,1,1),y=c(1,1,0), gp=gpar(col="black", fill="white"))#upper right triangle
					  grid.text("NA",x=0.8, y=0.6, gp=gpar(cex=0.7))
					  }
				  grid.polygon(x=c(0,0,1),y=c(1,0,0), gp=gpar(col="black", fill=i)) #bottom left triangle
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
