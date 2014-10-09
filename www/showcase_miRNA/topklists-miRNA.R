
## ----eval=FALSE----------------------------------------------------------
##     install.packages("TopKLists")


## ------------------------------------------------------------------------
    library(TopKLists)


## ------------------------------------------------------------------------
    res_51501 = read.table("data/GSE51501_results.txt", header=TRUE, quote="\"", stringsAsFactors=FALSE) 
    res_51504 = read.table("data/GSE51504_results.txt", header=TRUE, quote="\"", stringsAsFactors=FALSE) 
    res_51507 = read.table("data/GSE51507_results.txt", header=TRUE, quote="\"", stringsAsFactors=FALSE)


## ------------------------------------------------------------------------
    res_51507 = res_51507[complete.cases(res_51507),]
    res_51501 = res_51501[complete.cases(res_51501),]
    res_51504 = res_51504[complete.cases(res_51504),]


## ------------------------------------------------------------------------
    common = intersect(intersect(res_51501$mirname, res_51507$mirname), res_51504$mirname)
    length(common)


## ------------------------------------------------------------------------
    res_51507_common = res_51507[res_51507$mirname %in% common,]
    res_51501_common = res_51501[res_51501$mirname %in% common,]
    res_51504_common = res_51504[res_51504$mirname %in% common,]
    data_common = data.frame(HTS = res_51507_common$mirname, BeadChip = res_51501_common$mirname, NanoString = res_51504_common$mirname, stringsAsFactors=FALSE)
    dim(data_common)
    head(data_common)


## ----message=FALSE-------------------------------------------------------
    TopKListsGUI(data_common)
	
## ------------------------------------------------------------------------
	set.seed(123)
	res = j0.multi(data_common, d=40, v=22)	
	l1 =  as.character(data_common[1:res$maxK,1])
	l2 =  as.character(data_common[1:res$maxK,2])
	l3 =  as.character(data_common[1:res$maxK,3])
	input=list(l1,l2,l3)
	common=unique(unlist(input))
	space=list(common,common,common)
	outCEMC.kendall = CEMC(input, space, dm = "k")
	outCEMC.spearman = CEMC(input, space, dm = "s")

## ------------------------------------------------------------------------
    print(sessionInfo())


