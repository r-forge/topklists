# Algorithm for j0 estimation based on the method of Hall and Schimek (2012)
# The algorithm consists of an ordered sequence of "test stages" s1, s2, ....
# Stage sk is associated with an integer Jsk, which when k is odd, is a potential lower bound to j0
# Idata - input data is a vector of 0's and 1's

calculateDataSet <- function(lists, L, d, v, threshold) {
  message("StartCalculate")
  compared.lists <- list() #contains all pairwise compared lists (structure for aggmap)
  info <- matrix(ncol = 0, nrow = 3) #contains information about list names
  rownames(info) <- c("listname", "original listname", "ref-list or trunc-list")
  grayshaded.lists <- list() #contains information which element in the list has to be gray-shaded
  grayshaded.genes <- c() #contains all gray-shaded object-symbols or names
  temp.sumtrunclists <- list() #contains the summarized truncated lists (number of lists = L)
  summary.table <- matrix(nrow = 0, ncol = (4 + L)) #contains the summary-table
  venn.values <- list() #contains the Venn-lists (for viewing in the Venn-diagram) and the Venn-table (a Venn-diagram in table form)
  current.genesymbol_temp = c()
  
  
  
                                        #first step: calculate res.j0.temp
  res.j0.temp <- j0.multi(lists, d, v)
  res.temp <- as.matrix(res.j0.temp$L)
                                        #print(res.temp)
  max.j0.est <- res.j0.temp$maxK

##### adjusts v - in case the estimated max.j0.est is larger than nrow(lists)-v 
                                        #v_adj = min(nrow(lists)-max.j0.est, v, na.rm=T)

  temp = tapply(as.numeric(res.temp[,4]), res.temp[,1], function(x) max(x, na.rm = TRUE))
  
  if(sum(is.na(res.temp[,4]))<nrow(res.temp))
    {
      if (sum(temp!="-Inf")>1)
        {
          if (sum(is.na(res.temp[,4]))>0)
            {
              res.temp2 = res.temp[-which(is.na(res.temp[,4])),]
            } else {res.temp2 = res.temp}
        }else{
          if (sum(is.na(res.temp[,4]))>0)
            {res.temp2 = t(as.matrix(res.temp[-which(is.na(res.temp[,4])),]))
           } else {res.temp2 = t(as.matrix(res.temp))}
          
        }
      
      
      list_t = list()
      for (i in unique(res.temp2[,1]))
        {
          res.temp.temp = res.temp2[res.temp2[,1]==i,]
          if(!is.matrix(res.temp.temp)){res.temp.temp = t(as.matrix(res.temp.temp))}
          list_t[[i]] = cbind(res.temp.temp[,2][order(res.temp.temp[,4], decreasing=T)], res.temp.temp[,4][order(res.temp.temp[,4],decreasing=T)])
        }
      block_order = names(sort(unlist(lapply(list_t, FUN=function(x) max(as.numeric(x[,2]), na.rm=T))), decreasing=T))
                                        #print(block_order)
      inblock_list_order = lapply(list_t, FUN=function(x) x[,1][order(as.numeric(x[,2]), decreasing=T)])

      if (length(block_order)>1)
        {
          inblock_list_order_final = list()
          inblock_list_order_final[[block_order[1]]] = inblock_list_order[[block_order[1]]]
          block_order_temp = c()
          
          for (i in c(2:length(block_order)))
            {
              block_order_temp = c(block_order_temp, block_order[[i-1]]);
              inblock_list_order_final[[block_order[i]]] = setdiff(inblock_list_order[[block_order[i]]], block_order_temp)
            }
        }else{inblock_list_order_final = inblock_list_order}
      
      
      block_order_final = names(inblock_list_order_final)[which(lapply(inblock_list_order_final, length)!=0)]
      
##### condition if maximal estimated j0 is NA, then warning is returned #####

      if (max.j0.est!="-Inf")
        {	
          
###### building plotflow#####
          
          current.referencelist <- 0
                                        #iterate over all blocks (a block is a reference list with the corresponding truncation lists)
          
          for (first.list.name in block_order_final) {
            
                                        # selecting block
            temp2 = list_t[[first.list.name]]
            rownames(temp2) = list_t[[first.list.name]][,1]

                                        #get the object-symbols of the current reference list
                                        #gene.names.to.plot <- lists[, first.list.name][1:(as.numeric(temp2[inblock_list_order_final[[first.list.name]][1],2]))]
            gene.names.to.plot <- lists[, first.list.name][1:max.j0.est]
            gene.names.to.plot <- as.vector(gene.names.to.plot)
            
            current.referencelist <- current.referencelist + 1
            compared.lists[[paste("R", current.referencelist, sep = "")]] <- gene.names.to.plot
                                        #		temp.sumtrunclists[[inblock_list_order_final[[first.list.name]][1]]] <- lists[, inblock_list_order_final[[first.list.name]][1]][1:(as.numeric(temp2[inblock_list_order_final[[first.list.name]][1],2]))]
                                        #		info <- cbind(info, c(paste("R", current.referencelist, sep = ""), inblock_list_order_final[[1]][1], "R"))
            temp.sumtrunclists[[first.list.name]] <- gene.names.to.plot
            info <- cbind(info, c(paste("R", current.referencelist, sep = ""), first.list.name, "R"))
            message(temp.sumtrunclists)

            current.trunclist <- 0
            grayshaded.lists[[paste("R", current.referencelist, sep = "")]] <- rep(FALSE, length(gene.names.to.plot))		
            temp.countgray <- matrix(ncol = 2, nrow = length(gene.names.to.plot), data = 0)
            
                                        #iterate over the truncated lists of the current block
                                        #		if (length(inblock_list_order_final[[1]])>1){
            for (l in c(1:length(inblock_list_order_final[[first.list.name]]))) {
              message(paste(l, "\n"))
              n.genes.to.plot <- as.numeric(temp2[inblock_list_order_final[[first.list.name]][l],2])
              
              temp.distances = abs(c(1:length(gene.names.to.plot)) - match(gene.names.to.plot, as.character(lists[, inblock_list_order_final[[first.list.name]][l]])))
              temp.countgray[,2] = temp.countgray[,2]+c(temp.distances<=d)

              temp.sumtrunclists[[inblock_list_order_final[[first.list.name]][l]]] <- lists[, inblock_list_order_final[[first.list.name]][l]][1:(as.numeric(temp2[inblock_list_order_final[[first.list.name]][l],2]))]

                                        #check for gray-shade of an element in the truncated list
                                        #temp.grayshade = match(gene.names.to.plot, as.character(lists[, inblock_list_order_final[[first.list.name]][l]]))<(n.genes.to.plot)
              temp.grayshade = temp.distances<=d
                                        #add the truncated list
              current.trunclist <- current.trunclist + 1
              compared.lists[[paste("R", current.referencelist, "_T", current.trunclist, sep = "")]] <- temp.distances
              grayshaded.lists[[paste("R", current.referencelist, "_T", current.trunclist, sep = "")]] <- temp.grayshade
              info <- cbind(info, c(paste("R", current.referencelist, "_T", current.trunclist, sep = ""), inblock_list_order_final[[first.list.name]][l], "T"))
              
            }# end for l
            
                                        #calculate if respective object-symbol of the reference list has to be gray-shaded
            temp.percentage = apply(as.matrix(temp.countgray[,-1]),1,sum)/c(length(inblock_list_order_final[[first.list.name]]))*100
            grayshaded.lists[[paste("R", current.referencelist, sep = "")]] = temp.percentage >= threshold

                                        #add object-symbol to a new list which contains all gray-shaded object-symbols (add only if it is not already in the list)
                                        #print(lists[, first.list.name][which(temp.percentage >= threshold)])
            grayshaded.genes = union(grayshaded.genes, lists[, first.list.name][which(temp.percentage >= threshold)])
          }# end for first.list.name
          


                                        #having all the necessary information, calculate the summary-table
          colnames(summary.table) <- c(names(lists), "Rank sum", "Object order", "Freq in input lists", "Freq in truncated lists")
          
          for (j in 1:length(grayshaded.genes)) {
            current.genesymbol <- grayshaded.genes[j]
            
                                        #get the positions of the current object-symbol in the input lists
            temp.positions <- rep(NA,L)
            for (q in 1:L) {
              temp.positions[q] <- match(current.genesymbol, lists[,names(lists)[q]])
            }
            
                                        #calculate the rank sum
            temp.ranksum <- 0
            temp.missingvalues <- 0
            for (q in 1:L) {
              if (is.na(temp.positions[q])) {
                temp.missingvalues <- temp.missingvalues + 1
              }
            }
                                        #if one or more rank positions are unavaliable (e.g. an object-symbol does not exist in a list), interpolate the rank sum
            if (temp.missingvalues > 0) {
              temp.meanrank <- 0
              temp.partialranksum <- 0
                                        #get the rank sum of the valid rank-positions
              for (q in 1:L) {
                if (!is.na(temp.positions[q])) {
                  temp.partialranksum <- temp.partialranksum + temp.positions[q]
                }
              }
                                        #calculate the mean rank-position for the unavailable rank-positions
              temp.meanrank <- round(temp.partialranksum / (L - temp.missingvalues))
              temp.ranksum <- temp.partialranksum + (temp.meanrank * temp.missingvalues)
            } else {
              temp.ranksum <- sum(temp.positions)
            }
            
                                        #calculate the frequency in the input lists
            temp.freqinput <- L - temp.missingvalues
            
                                        #calculate the frequency in the summarized truncated lists
            temp.freqtrunc <- 0	
            for (curr.listname in names(temp.sumtrunclists)) {
              if (current.genesymbol %in% temp.sumtrunclists[[curr.listname]]) {
                temp.freqtrunc <- temp.freqtrunc + 1
              }
            }		
            
                                        #add the calculated row (of the current object) to the summary-table
            current.genesymbol_temp = c(current.genesymbol_temp, current.genesymbol)
            summary.table <- rbind(summary.table, c(temp.positions, temp.ranksum, NA, temp.freqinput, temp.freqtrunc))
          }# end for j
          
                                        #conversion of the summary table into a data frame so that the rankings are given as numbers not as char, otherwise the ordering is wrong
          
          summary.table.final = data.frame(Object=current.genesymbol_temp,summary.table, stringsAsFactors=FALSE)
          
          
                                        #the last task for creation of the summary-table is to order the object-symbols according to their rank-sum
          temp.counter <- 0
          for (curr.genepos in order(summary.table.final[,(2 + L)])) {
            temp.counter <- temp.counter + 1
            summary.table.final[temp.counter,(3+L)] <- summary.table.final[curr.genepos,1]
          }
          
                                        #calculate the Venn-lists (to view the Venn-diagram) and the Venn-table
                                        #the calculation takes place only for L between 2 and 4 (a Venn-diagram for L > 4 cannot be properly arranged)
          venn.values <- .calculateVennValues(summary.table.final[,1:(L + 1)], L)
          
                                        #combine all necessary objects into one single list
          truncated.lists <- list()
          truncated.lists$comparedLists <- compared.lists
          truncated.lists$info <- info
          truncated.lists$grayshadedLists <- grayshaded.lists
          truncated.lists$summarytable <- summary.table.final
          truncated.lists$vennlists <- venn.values$vennlists
          truncated.lists$venntable <- venn.values$venntable
          truncated.lists$v <- v
          truncated.lists$Ntoplot<-sum(unlist(lapply(inblock_list_order_final,length)))+sum(unlist(lapply(inblock_list_order_final,length))>0)
          truncated.lists$Idata <- res.j0.temp$Idata
          return(truncated.lists)
          #message(truncated.lists)
        }
    } else {
      message(paste("!!!...For selected delta, the top L list cannot be estimated (little or no overlap)!!!", "\n"))
      return(truncated.lists=NULL)
    } # end if if (max.j0.est)
}# end of function calculateDataSet

compute.stream<-function(Idata, const=0.251, v, r=1.2)
{
if(sum(Idata, na.rm=T)==length(Idata)) 
{
return(list(j0_est=length( ), reason.break="Idata is identity", Js=NA, v=v))

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

#Method of Hall and Schimek (2012) adapted to the situation of multiple ranked lists 
#Lists have the same length and comprise no missing values
#Inputs:
#Data - input matrix, each column represents one list
#delta - the maximal distance between the rank positions of an object in a pair of lists
#v - value of tuning parameter nu used for j0 estimation
#Outputs: maximal estimated j0 based on all combinations of 2 lists

#Data = lists
#delta=10
#v=10

j0.multi<-function(lists,d,v) {
  maxK=0
  L = c()
  Idata_ID = c()
  names_idata = c()
  for (i in 1:ncol(lists)){
  for (j in 1:ncol(lists)){
    if (i!=j) {
	ID = prepareIdata(lists[,c(i,j)],d=d)
        Idata_ID=cbind(Idata_ID, ID$Idata)
	names_idata = c(names_idata, paste("L",i,"_","L",j,sep=""))
	J = compute.stream(ID$Idata,v=v)$j0_est
	L = rbind(L, cbind(names(lists)[i], names(lists)[j], v,J,d))
	}# end for if
      }# end for j
    }# end for i
colnames(Idata_ID) = names_idata
L = data.frame(L, stringsAsFactors=F)
names(L) = c("list1", "list2", "v", "j0_est","delta")
if(sum(is.na(L$j0_est))<nrow(L)){maxK = max(as.numeric(L$j0_est), na.rm=T)}else{maxK=NA}
return(list(maxK=maxK,L=L, Idata=Idata_ID))
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

# Function to prepare Idata from several rankings of several assessors
# x - data matrix, where columns represent the rank order of objects from two different assessors and the rows represent object names (symbols)
# delta - the maximal distance between the rank positions of an object in a pair of lists
# num.omit - the maximal number of ommited objects from the analysis 
#
# The result is an object of type "Idata", which is a list containing Idata and the information on the distance delta.

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
