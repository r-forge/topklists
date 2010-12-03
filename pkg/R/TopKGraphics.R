
## Heatmap solution to HS rank order graphics representation
## By Eva Budinska, 3.12.2010
## Input variables 
#K- number of lists
# N maximal number of genes
# delta=10 maximal distance between genes used for j0 estimation
# v - nu value used for j0 estimation
# res.temp - data frame, result of HS algorithm for all pairs of lists, 4 columns, first and second column contain lists names that were compared, third column selected nu and fourth column estimated j0
# lists data frame containing three orderesd lists of genes that were compared (the input to HS algorithm) - column names obligatory

aggmap<-function(lists, res.j0.temp, N, K, delta, v)
{

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
library(grid)
grid.newpage()

# defining basic layout with two rows and one column
pushViewport(viewport(layout=grid.layout(2, 1, widths=c(1,1),heights=c(0.1,0.9))))

##first row - defining header (delta, nu and color key level) 
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
