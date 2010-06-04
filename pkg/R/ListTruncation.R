# method for j0 estimation
# the algorithm consists of an ordered sequence of "test stages" s1, s2....
# stage sk is associated with an integer Jsk, which when k is odd, is a potential lower bound to j0
# Idata - input data is a vector of 0's and 1's

compute.stream <- function(Idata, const=0.6, v, r=1.2){
#library(sma)
#source("pjplus.R")
#source("pjminus.R")
#source("moderate.deviation.R")

zv 		= .moderate.deviation(const, v)
Js 		= c()
k 		= 1
pj.plus 	= c()
pj.minus 	= c()
h		= 0

###########
#if (sum(Idata)==length(Idata)){j0_est=length(Idata)}
reason.break="NA"
repeat
{
if (is.odd(k))
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
		
			if ((j+v[k]-1) >= length(Idata)) {reason.break="(j+v[k]-1) >= length(Idata)"
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
		print(Js)

		}# end if

	} else{
		h = Js[k-1]-trunc(r*v.last)+1;
		j = h-1;
		#print(paste("j=",j, sep=""))
		if (j< length(Idata))
		{
			repeat 
			{
			j = j+1
		
			if ((j+v.last-1) >= length(Idata)) {reason.break="(j+v.last-1) >= length(Idata)"
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
		print(Js)
		# breaking the repeat loop condition 2
#		if ((k-3)>=1) {if (Js[k-2]==Js[k] & Js[k-1]==Js[k-3]) {reason.break="Js[k-2]==Js[k] & Js[k-1]==Js[k-3]";
#							break} }
		if ((k-3)>=1) 
			{
			if (Js[k-2]==Js[k] & Js[k-1]==Js[k-3]) 
				{
				if(abs(Js[k-1]-Js[k])<=2)
					{reason.break="Js[k-2]==Js[k] & Js[k-1]==Js[k-3]";
#					break;
				}else{
					j=j/2
					v.last = ceiling(v.last/2)}
				} 
			
			}


		}# end if
	}#end else
	if(reason.break!="NA"){break;}
} # end is.odd

if (is.even(k))
{
	h = Js[k-1]+trunc(r*v.last)
	j = h

	if (j>1)
	{
		repeat
		{
		j = j-1
		print(j)
		print(reason.break)

		# computing pjminus
		if ((j-v.last+1) <= 0){reason.break="(j-v.last+1) <= 0"
					break}
		
		pj.minus 	= .pjminus(Idata, v.last, j)
		#print(pj.minus)
		#print(Idata[(j-v.last+1):j])
			
		# testing
		if ((pj.minus-0.5)>zv) {#print("pj.minus-0.5>zv")
						break}	
		}# end repeat

	Js[k] = j
	print(Js)
	if(reason.break!="NA"){break;}
	# breaking the repeat loop condition 3
	#if (Js[k]-trunc(r*v)+1<=1) {reason.break="Js[k]-r*v<=1"; break}
	 if (Js[k]-trunc(r*v.last)+1<=1) {v.temp=v.last; 
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
print(paste("k=",k, sep=""))

}#end repeat

if(v.last<=1 | reason.break=="(j-v.last+1) <= 0"){j0_est=0} else {if (is.odd(k) & k>2) {j0_est = ceiling(Js[k-2]+0.5*v[k-2]-1)}else if (is.even(k) & k>1){j0_est = ceiling(Js[k-1]+0.5*v[k-1]-1)} else{j0_est=length(Idata)}}

return(list(j0_est=j0_est, reason.break=reason.break, Js=Js, v.vector=v))
}

#Method of M. G. Schimek extended to count K(j0_est) for each pair of several ranked lists. 
#Lists have same length and have no missing values.
#Inputs:
#Data - input matrix, each coloum presents one list
#m.dist - the maximal distance of the objects in ranking
#v - parameter for estimating j0
#outputs: maximal estimated J0 through all combinations of 2 lists
j0.multi<-function(Data,m.dist,v) {
  maxK=0
  K=matrix(NA,ncol(Data), ncol(Data))

  for (A in 1:ncol(Data)){
  for (B in 1:ncol(Data)){
    if (A!=B) {
	ID=prepareIdata(Data[,c(A,B)],type="second",m.dist=m.dist)
	#orig (wrong return for matrix)
	#J=compute.stream(ID$Idata,v=v)
	# modified by Karl Kugler
	#TODO: should output be j0.est
	J=compute.stream(ID$Idata,v=v)$j0_est
	K[A,B]=J
	if (J>maxK){
	  maxK<-J
	}
      }
    }
  }
 return(list(maxK,K))
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
# x 	- data matrix, where 
# # if type = first: first column represents rank order of genes from the first assessor, second column represents assignment to the category (0, 1; provided by second assessor)
# # if type = second:  columns represents rank order of genes from two different assessors, rows represents genes (clones)

# type - "first" or "second", represents type of the mechanism to create Idata
# genelist - a vector containing clone names. Must be the same length as the number of rows in x. 
## 		Order of genes in x must be the same as in genelist. 
# expr - whether the data are expression values (will be needed in further version of the algorithm) 
# m.dist - the maximal distance of the objects in ranking between assessors
# num.omit - the maximal number of ommited genes from the analysis 
#
# the result is a object of type "Idata", which is a list containing Idata and the information
# about the distance (m.dist).

prepareIdata <- function(x, type, genelist, expr = "FALSE", m.dist)
{
Idata = rep(0, nrow(x))
if (type == "second")
	{
	rank.diff = c()
	for (i in 1:nrow(x))
		{
		rank.diff = c(rank.diff, abs(i-which(x[,2]==x[i,1])))
		Idata[which(rank.diff<=m.dist)] = 1
		}
	}
if (type == "first")
	{
	Idata = x[2];
	}
return(list(Idata = Idata, m.dist = m.dist))
}

is.odd <- function(k){
  return(k %% 2 != 0)
}

is.even <- function(k){
  return(!is.odd(k))
}
