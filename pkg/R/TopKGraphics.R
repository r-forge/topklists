TopKListsGUI <- function(lists, autorange.delta = FALSE, override.errors = FALSE, aggmap.size = c(870, 440), aggmap.res = 100, venndiag.size = c(380, 420), venndiag.res = 70, gui.size = c(900, 810), directory = "TopKListsGUI_temp") {
# 	require(RGtk2)
# 	require(gWidgets)
# 	require(gWidgetsRGtk2)
# 	library(gplots)
	options("guiToolkit"="RGtk2")
	
	#delta_symbol = substitute(delta)
	#nu_symbol = expression(nu)

	#check if given lists are of type 'data.frame' and contains data
	if (!is.data.frame(lists)) {
		gmessage("ERROR:\nGiven list-object is not of type data.frame!\nProgram will exit.", title = "Incorrect input", icon = "error")
		return("TopKLists GUI Error: Given list-object is not of type data.frame!")
	} else if ((dim(lists)[1] < 2) | (dim(lists)[2] < 2)) {
		gmessage("ERROR:\nGiven list-object has eighter only one feature or only one list!\nProgram will exit.", title = "Incorrect input", icon = "error")
		return("TopKLists GUI Error: Given list-object has eighter only one feature or only one list!")
	}

	win <- gwindow("TopKLists GUI", visible = FALSE, width = gui.size[1], height = gui.size[2]) #main window
	maingroup <- ggroup(horizontal = FALSE, container = win, expand = TRUE) #main container
	wid <- environment() #contains all the widgets
	current.arguments <- environment() #contains all the arguments (N,K,v etc) of the current calculation

	#create input-area for arguments
	calcframe <- gframe("ARGUMENTS", pos = 0.5, container = maingroup)
	calcgroup <- ggroup(horizontal = FALSE, container = calcframe, expand = TRUE)
	agroup <- ggroup(horizontal = TRUE, container = calcgroup, expand = TRUE)
	a1frame <- gframe("Calculated from given list", container = agroup, expand = TRUE)
	a2frame <- gframe("Choose range for delta (depends on v)", container = agroup, expand = TRUE)
	a3frame <- gframe("Threshold", container = agroup, expand = TRUE)

	wid$calculate <- gbutton("Calculate", container = calcgroup, handler = function(h,...) { calculate.all.truncationlists() })	
	svalue(calcgroup) <- 10

	#layout of the first argument-input-area
	tbl1 <- glayout(spacing = 10,container = a1frame)
	tbl1[2,2, anchor = c(1,0)] <- glabel("N", container = tbl1)
	tbl1[2,3] <- wid$N <- gspinbutton(from = 0, to = dim(lists)[1], by = 1, value = dim(lists)[1], container = tbl1)
	tbl1[2,4, anchor = c(-1,0)] <- glabel("maximal number of features" ,container = tbl1)
	tbl1[3,2, anchor = c(1,0)] <- glabel("L", container = tbl1)
	tbl1[3,3] <- wid$K <- gspinbutton(from = 0, to = dim(lists)[2], by = 1, value = dim(lists)[2], container = tbl1)
	tbl1[3,4, anchor = c(-1,0)] <- glabel("number of lists", container = tbl1)

	#layout of the second argument-input-area
	tbl2 <- glayout(spacing = 10, container = a2frame)
	tbl2[2,2, anchor = c(1,0)] <- glabel(expression(nu), container = tbl2)
	tbl2[2,3] <- wid$v <- gspinbutton(from = 0, to = dim(lists)[1], by = 1, value = 0, container = tbl2, handler = function(h,...) { update.deltarange() })
	tbl2[2,4:5, anchor = c(1,0)] <- glabel("value used for j0 estimation", container = tbl2)
	tbl2[3,2, anchor = c(1,0)] = glabel(expression(delta), container = tbl2)
	tbl2[3,3] <- wid$delta.start <- gspinbutton(from = 0, to = dim(lists)[1], by = 1, value = 0, container = tbl2)
	if (autorange.delta) {
		enabled(wid$delta.start) <- FALSE
	}
	tbl2[3,4, anchor = c(0,0)] = glabel("<- between ->", container = tbl2)
	tbl2[3,5] <- wid$delta.stop <- gspinbutton(from = 0, to = dim(lists)[1], by = 1, value = 0, container = tbl2)
	if (autorange.delta) {
		enabled(wid$delta.stop) <- FALSE
	}

	#layout of the third argument-input-area
	tbl3 <- glayout(spacing = 10, container = a3frame)
	tbl3[2,2, anchor = c(1,0)] <- glabel("min", container = tbl3)
	tbl3[2,3] <- wid$threshold <- gedit(text = "50", width = 3, container = tbl3)
	tbl3[2,4, anchor = c(-1,0)] <- glabel("%" ,container = tbl3)
	tbl3[3,2:4, anchor = c(1,0)] <- glabel("to gray-shade gene symbol", container = tbl3)

	#create the delta-slider for choosing the desired delta-value (and to show the calculated data for the selected delta-value)
	sldrframe <- gframe("DELTA-SLIDER", pos = 0.5, container = maingroup)
	#the parameters 'from' and 'to' are dummy-values because the slider will be updated after the calculation
	wid$delta.slider <- gslider(from = 1, to = 2, by = 1, container = sldrframe, expand = TRUE, handler = function(h,...) { load.data() })
	enabled(wid$delta.slider) <- FALSE

	#create three three tabs for presenting the results of the calculation
	nb <- gnotebook(container = maingroup, expand = TRUE)
	aggmapg <- ggroup(horizontal = FALSE, container = nb, label = "Aggregation map")
	summtblg <- ggroup(horizontal = FALSE, container = nb, label = "Summary table")
	venng <- ggroup(horizontal = TRUE, container = nb, label = "Venn-diagram & Venn-table")
	svalue(nb) <- 1 #set the first tab as the selected tab

	#create a status-bar to inform the user of current events
	wid$status <- gtkStatusbar()
	add(maingroup, wid$status)
	info <- wid$status$getContextId("info")
	wid$status$push(info, "GUI started")

	#create a progress-bar to show the progress of the calculation
	wid$progress <- gtkProgressBar(show = TRUE)
	add(maingroup, wid$progress)	

	#update the allowed range for delta (based on the v-value)
	update.deltarange <- function() {
		#check if allowed range for delta should be updated
		if (autorange.delta) {
			deltarange <- .determineDeltaRange(lists, as.numeric(svalue(wid$N)), as.numeric(svalue(wid$K)), as.numeric(svalue(wid$v)))
			svalue(wid$delta.start) <- deltarange[1]
			svalue(wid$delta.stop) <- deltarange[2]
			enabled(wid$delta.start) <- TRUE
			enabled(wid$delta.stop) <- TRUE
			wid$status$push(info, expression(paste("Recalculated allowed range for delta with ",nu," =", ss), list(ss=svalue(wid$v))))
		}
	}

	#update the calculated data set (aggmap, summary-table and venn) in the three tabs
	load.data <- function() {		
		#update is only performed when delta-slider is enabled
		if (enabled(wid$delta.slider)) {
			#check if calculated data set can be shown in the GUI (no file = error in the calculation for this delta-value)
			if (file.exists(paste(directory, "/N" , current.arguments$val.N, "_L", current.arguments$val.K, "_delta", svalue(wid$delta.slider), "_v", current.arguments$val.v, "_thrshld", current.arguments$val.threshold, ".Rdata", sep = ""))) {
				#load the calculated data set
				load(paste(directory, "/N", current.arguments$val.N, "_L", current.arguments$val.K, "_delta", svalue(wid$delta.slider), "_v", current.arguments$val.v, "_thrshld", current.arguments$val.threshold, ".Rdata", sep = ""))

				#first tab: show aggmap
				if (!is.null(wid$aggmap)) {
					#delete current aggmap-image
					delete(aggmapg, wid$aggmap)
				}
				wid$aggmap <- gimage(paste(directory, "/aggmap_N", current.arguments$val.N, "_L", current.arguments$val.K, "_delta", svalue(wid$delta.slider), "_v", current.arguments$val.v, "_thrshld", current.arguments$val.threshold, ".png", sep = ""), container = aggmapg)

				#second tab: show summary-table
				if (!is.null(wid$sumtable)) {
					#delete current summary-table
					delete(summtblg, wid$sumtable)
				}
				wid$sumtable <- gtable(truncated.lists$summarytable, container = summtblg, expand = TRUE)

				#third tab: show venn-diagram and venn-table (only if K <= 4)
				if (current.arguments$val.K <= 4) {			
					if (!is.null(wid$venndiag) & !is.null(wid$venntbl)) {
						#delete current venn-diagram and venn-table
						delete(venng, wid$venndiag)
						delete(venng, wid$venntbl)
					}
					if (!is.null(wid$nodraw) & !is.null(wid$nodrawimage)) {
						#delete the icon and message
						delete(venng, wid$nodrawimage)
						delete(venng, wid$nodraw)
					}
					wid$venndiag <- gimage(paste(directory, "/venn_N", current.arguments$val.N, "_L", current.arguments$val.K, "_delta", svalue(wid$delta.slider), "_v", current.arguments$val.v, "_thrshld", current.arguments$val.threshold, ".png", sep = ""), container = venng)
					wid$venntbl <- gtable(truncated.lists$venntable, container = venng, expand = TRUE)
				} else {
					#output message that venn-diagram and venn-table are not shown for K > 4
					if (is.null(wid$nodrawimage) & is.null(wid$nodraw)) {
						wid$nodrawimage <- gimage("info", dirname = "stock", container = venng)
						wid$nodraw <- glabel("There are more than four input lists.\nVenn-diagram and Venn-table are only drawn to a maximum of four input lists.", container = venng)
					}
				}
			} else {
				gmessage("There is no data to show in the GUI!\nThis is caused by an error in the calculation.\nYou may have to adjust the arguments!", title = "Loading data failed", icon = "error")
			}
		}
	}
	
	#calculates a data set (aggmap, summary-table and venn) for each delta in the given range for delta
	calculate.all.truncationlists <- function() {
		#check if entered values are valid
		if ((as.numeric(svalue(wid$delta.start)) < 0) || (as.numeric(svalue(wid$delta.stop)) < 0) || (svalue(wid$delta.start) > svalue(wid$delta.stop))) {
			gmessage("Invalid range for delta entered.\nPlease check your input.", title = "Invalid range", icon = "error")
			return()
		}
		temp.threshold <- as.numeric(svalue(wid$threshold))
		if (is.na(temp.threshold) || (temp.threshold < 0) || (temp.threshold > 100)) {
			gmessage("The entered threshold is not a valid number!\nPlease enter a threshold between 0 and 100.", title = "Threshold invalid", icon = "error")
			return()
		}

		wid$status$push(info, "Begin calculation of data sets for range of delta")
		gtkProgressBarSetFraction(wid$progress, 0) #set progress-bar to start-value
		gtkMainIterationDo(FALSE)

		#disable delta-slider to prevent sliding while calculating
		enabled(wid$delta.slider) <- FALSE

		#delete old files in destination directory
		unlink(directory, recursive = TRUE)
		dir.create(directory)	

		#calculate data set for each delta in the given range for delta
		temp.tocalc <- (as.numeric(svalue(wid$delta.stop)) - as.numeric(svalue(wid$delta.start))) + 1
		temp.loopcounter <- 0
		wid$error <- FALSE
		for (i in as.numeric(svalue(wid$delta.start)):as.numeric(svalue(wid$delta.stop))) {
			gtkMainIterationDo(FALSE)
			#try to calculate data set for current delta
			tryCatch({
				truncated.lists <- .calculateDataSet(lists, as.numeric(svalue(wid$K)), i, as.numeric(svalue(wid$v)), as.numeric(svalue(wid$threshold)))
				#save calculated truncated lists for current delta in the destination directory
				save(truncated.lists, file = paste(directory, "/N", as.numeric(svalue(wid$N)), "_L", as.numeric(svalue(wid$K)), "_delta", i, "_v", as.numeric(svalue(wid$v)), "_thrshld", as.numeric(svalue(wid$threshold)), ".Rdata", sep = ""))
				#draw the aggmap and save it as png-image in the destination directory
				png(filename = paste(directory, "/aggmap_N", as.numeric(svalue(wid$N)), "_L", as.numeric(svalue(wid$K)), "_delta", i, "_v", as.numeric(svalue(wid$v)), "_thrshld", as.numeric(svalue(wid$threshold)), ".png", sep = ""), width = aggmap.size[1], height = aggmap.size[2], res = aggmap.res)
				.drawAggmap(truncated.lists, N=as.numeric(svalue(wid$N)), K=as.numeric(svalue(wid$K)), d=i, v=as.numeric(svalue(wid$v)), threshold=as.numeric(svalue(wid$threshold)), lists)
				dev.off()
				#draw the venn-diagram and save it as png-image in the destination directory (only if K is between 2 and 4) & if there is a j0 estimated
				if((as.numeric(svalue(wid$K)) >= 2) & (as.numeric(svalue(wid$K)) <= 4) & !is.null(truncated.lists)) {
					png(filename = paste(directory, "/venn_N", as.numeric(svalue(wid$N)), "_L", as.numeric(svalue(wid$K)), "_delta", i, "_v", as.numeric(svalue(wid$v)), "_thrshld", as.numeric(svalue(wid$threshold)), ".png", sep = ""), width = venndiag.size[1], height = venndiag.size[2], bg = "transparent", res = venndiag.res)
					venn(truncated.lists$vennlists)
					dev.off()
				}else{
				png(filename = paste(directory, "/venn_N", as.numeric(svalue(wid$N)), "_L", as.numeric(svalue(wid$K)), "_delta", i, "_v", as.numeric(svalue(wid$v)), "_thrshld", as.numeric(svalue(wid$threshold)), ".png", sep = ""), width = venndiag.size[1], height = venndiag.size[2], bg = "transparent", res = venndiag.res)
				plot(1,1, type="n", axes=FALSE)
				text(1,1,"No overlap on selected parameters")
				dev.off()

				}
			}, error = function(e) {
				#catch errors while calculating the data set
				wid$status$push(info, paste("Error in the current calculation for delta =", i))
				wid$error <- TRUE
				if(!override.errors) {
					return()
				}
			})# end tryCatch
			if (wid$error & !override.errors) {
				#exit the loop
				break()
			}

			#update progress-bar
			temp.loopcounter <- temp.loopcounter + 1
			temp.fraction <- (1 / temp.tocalc) * temp.loopcounter
			gtkProgressBarSetFraction(wid$progress, temp.fraction)
			gtkMainIterationDo(FALSE)
		}

		#if errors are overridden, the calculation continues (but there is no data to view in the GUI for this delta value)
		if (wid$error) {
			if(override.errors) {
				wid$status$push(info, "Finished calculation of data sets for range of delta (with errors occurred)")
			} else {
				enabled(wid$delta.slider) <- FALSE
				return()
			}
			wid$error <- FALSE
		} else {		
			wid$status$push(info, "Finished calculation of data sets for range of delta")
		}

		#update the range of the delta-slider
		wid$delta.slider[] <- seq(as.numeric(svalue(wid$delta.start)), as.numeric(svalue(wid$delta.stop)), by = 1)

		#set the start position of the slider and enable slider
		if (temp.tocalc > 2) {
			svalue(wid$delta.slider) <- as.numeric(svalue(wid$delta.start)) + round((temp.tocalc - 1) / 2)
		} else {
			svalue(wid$delta.slider) <- as.numeric(svalue(wid$delta.start))
		}
		enabled(wid$delta.slider) <- TRUE

		#save the current arguments
		#(this prevents errors when sliding through the delta values - otherwise changing the arguments in the GUI could cause errors when showing the calculated data)
		current.arguments$val.N <- as.numeric(svalue(wid$N))
		current.arguments$val.K <- as.numeric(svalue(wid$K))
		current.arguments$val.v <- as.numeric(svalue(wid$v))
		current.arguments$val.threshold <- as.numeric(svalue(wid$threshold))
		
		#show the data in the tabs
		load.data()
	}#end calculate.all.truncationlists
	
	#show the GUI
	visible(win) <- TRUE
}



.calculateDataSet <- function(lists, K, d, v, threshold) {
	print("StartCalculate")
	compared.lists <- list() #contains all pairwise compared lists (structure for aggmap)
	info <- matrix(ncol = 0, nrow = 3) #contains information about list names
	rownames(info) <- c("listname", "original listname", "ref-list or trunc-list")
	grayshaded.lists <- list() #contains information which element in the list has to be gray-shaded
	grayshaded.genes <- c() #contains all gray-shaded gene-symbols
	temp.sumtrunclists <- list() #contains the summarized truncated lists (number of lists = K)
	summary.table <- matrix(nrow = 0, ncol = (4 + K)) #contains the summary-table
	venn.values <- list() #contains the venn-lists (to view in the venn-diagram) and the venn-table (a venn-diagram as a table)
	current.genesymbol_temp = c()
	
	
	
	#first step: calculate res.j0.temp
	res.j0.temp <- j0.multi(lists, d, v)
	res.temp <- as.matrix(res.j0.temp$K)
	#print(res.temp)
	max.j0.est <- res.j0.temp$maxK

	##### adjusting v - in case the estimated max.j0.est is larger than nrow(lists)-v #####EBcorr
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
	
		##### adding condition if maximal estimated j0 is NA, then return warning #####

	   if (max.j0.est!="-Inf")#####EBcorr
	   {	
	
		###### building plotflow#####
		
		current.referencelist <- 0
		#iterate over all blocks (a block is a reference list with the corresponding truncation lists)
	
		for (first.list.name in block_order_final) {
 
	 	# selecting block
	 	temp2 = list_t[[first.list.name]]
	 	rownames(temp2) = list_t[[first.list.name]][,1]

		#get the gene-symbols of the current reference list
		#gene.names.to.plot <- lists[, first.list.name][1:(as.numeric(temp2[inblock_list_order_final[[first.list.name]][1],2]))]
		gene.names.to.plot <- lists[, first.list.name][1:max.j0.est]
		gene.names.to.plot <- as.vector(gene.names.to.plot)
		
		current.referencelist <- current.referencelist + 1
		compared.lists[[paste("R", current.referencelist, sep = "")]] <- gene.names.to.plot
#		temp.sumtrunclists[[inblock_list_order_final[[first.list.name]][1]]] <- lists[, inblock_list_order_final[[first.list.name]][1]][1:(as.numeric(temp2[inblock_list_order_final[[first.list.name]][1],2]))]
#		info <- cbind(info, c(paste("R", current.referencelist, sep = ""), inblock_list_order_final[[1]][1], "R"))
		temp.sumtrunclists[[first.list.name]] <- gene.names.to.plot
		info <- cbind(info, c(paste("R", current.referencelist, sep = ""), first.list.name, "R"))
		print(temp.sumtrunclists)

		current.trunclist <- 0
		grayshaded.lists[[paste("R", current.referencelist, sep = "")]] <- rep(FALSE, length(gene.names.to.plot))		
		temp.countgray <- matrix(ncol = 2, nrow = length(gene.names.to.plot), data = 0)
		
		#iterate over the truncated lists of the current block
#		if (length(inblock_list_order_final[[1]])>1){
		for (k in c(1:length(inblock_list_order_final[[first.list.name]]))) {
		cat(paste(k, "\n"))
			n.genes.to.plot <- as.numeric(temp2[inblock_list_order_final[[first.list.name]][k],2])
			
			temp.distances = abs(c(1:length(gene.names.to.plot)) - match(gene.names.to.plot, as.character(lists[, inblock_list_order_final[[first.list.name]][k]])))
			temp.countgray[,2] = temp.countgray[,2]+c(temp.distances<=d)

			temp.sumtrunclists[[inblock_list_order_final[[first.list.name]][k]]] <- lists[, inblock_list_order_final[[first.list.name]][k]][1:(as.numeric(temp2[inblock_list_order_final[[first.list.name]][k],2]))]

#			#check if to gray-shade the element in the truncated list
			#temp.grayshade = match(gene.names.to.plot, as.character(lists[, inblock_list_order_final[[first.list.name]][k]]))<(n.genes.to.plot)
			temp.grayshade = temp.distances<=d
			#add the truncated list
			current.trunclist <- current.trunclist + 1
			compared.lists[[paste("R", current.referencelist, "_T", current.trunclist, sep = "")]] <- temp.distances
			grayshaded.lists[[paste("R", current.referencelist, "_T", current.trunclist, sep = "")]] <- temp.grayshade
			info <- cbind(info, c(paste("R", current.referencelist, "_T", current.trunclist, sep = ""), inblock_list_order_final[[first.list.name]][k], "T"))
	
			}# end for k
	
			#calculate if respective gene-symbol of the reference list has to be gray-shaded
			temp.percentage = apply(as.matrix(temp.countgray[,-1]),1,sum)/c(length(inblock_list_order_final[[first.list.name]]))*100
			grayshaded.lists[[paste("R", current.referencelist, sep = "")]] = temp.percentage >= threshold

			#add gene-symbol to a new list which contains all gray-shaded gene-symbols (add only if it is not already in the list)
			#print(lists[, first.list.name][which(temp.percentage >= threshold)])
			grayshaded.genes = union(grayshaded.genes, lists[, first.list.name][which(temp.percentage >= threshold)])
		}# end for first.list.name
	


			#having all the necessary information, calculate the summary-table
			colnames(summary.table) <- c(names(lists), "Rank sum", "Object order", "Freq in input lists", "Freq in truncated lists")
		
		for (j in 1:length(grayshaded.genes)) {
		current.genesymbol <- grayshaded.genes[j]
		
		#get the positions of the current gene-symbol in the input lists
			temp.positions <- rep(NA,K)
			for (q in 1:K) {
				temp.positions[q] <- match(current.genesymbol, lists[,names(lists)[q]])
			}
		
			#calculate the rank sum
			temp.ranksum <- 0
			temp.missingvalues <- 0
			for (q in 1:K) {
				if (is.na(temp.positions[q])) {
					temp.missingvalues <- temp.missingvalues + 1
				}
			}
			#if one or more rank-positions are not avaliable (e.g. a gene-symbol does not exist in a list), interpolate the rank sum
			if (temp.missingvalues > 0) {
				temp.meanrank <- 0
				temp.partialranksum <- 0
				#get the rank sum of the valid rank-positions
				for (q in 1:K) {
					if (!is.na(temp.positions[q])) {
						temp.partialranksum <- temp.partialranksum + temp.positions[q]
					}
				}
				#calculate the mean rank-position for the not available rank-positions
				temp.meanrank <- round(temp.partialranksum / (K - temp.missingvalues))
				temp.ranksum <- temp.partialranksum + (temp.meanrank * temp.missingvalues)
			} else {
				temp.ranksum <- sum(temp.positions)
			}
		
			#calculate the frequency in the input lists
			temp.freqinput <- K - temp.missingvalues
		
			#calculate the frequency in the summarized truncated lists
			temp.freqtrunc <- 0	
			for (curr.listname in names(temp.sumtrunclists)) {
				if (current.genesymbol %in% temp.sumtrunclists[[curr.listname]]) {
					temp.freqtrunc <- temp.freqtrunc + 1
				}
			}		
		
			#add the calculated row (of current object) to the summary-table
			current.genesymbol_temp = c(current.genesymbol_temp, current.genesymbol)
			summary.table <- rbind(summary.table, c(temp.positions, temp.ranksum, NA, temp.freqinput, temp.freqtrunc))
			}# end for j
		
		#converting the summary table into data frame so that the rankings are as numbers not as char, otherwise the ordering is wrong
		
		summary.table.final = data.frame(Object=current.genesymbol_temp,summary.table, stringsAsFactors=FALSE)
				
		
		#the last task for creating the summary-table is to order the gene-symbols according to their rank-sum
		temp.counter <- 0
		for (curr.genepos in order(summary.table.final[,(2 + K)])) {
			temp.counter <- temp.counter + 1
			summary.table.final[temp.counter,(3+K)] <- summary.table.final[curr.genepos,1]
		}
	
		#calculate the venn-lists (to view the venn-diagram) and the venn-table
		#the calculation takes place only fo K between 2 and 4 (a venn-diagram for K > 4 is not clearly arranged)
		venn.values <- calculateVennValues(summary.table.final[,1:(K + 1)], K)
	
		#combine all necessary objects into one single list
		truncated.lists <- list()
		truncated.lists$comparedLists <- compared.lists
		truncated.lists$info <- info
		truncated.lists$grayshadedLists <- grayshaded.lists
		truncated.lists$summarytable <- summary.table.final
		truncated.lists$vennlists <- venn.values$vennlists
		truncated.lists$venntable <- venn.values$venntable
		truncated.lists$v <- v####EBcorr
		truncated.lists$Ntoplot<-sum(unlist(lapply(inblock_list_order_final,length)))+sum(unlist(lapply(inblock_list_order_final,length))>0)####EBcorr
	
		return(truncated.lists)
        	print(truncated.lists)
		}
	   } else {
		cat(paste("!!!...For selected delta, the top K list cannot be estimated (little or no overlap)!!!", "\n"))
		return(truncated.lists=NULL)
} # end if if (max.j0.est)	#####EBcorr
}

.calculateVennValues <- function(genetable, K) {
	#check if K is between 2 and 4 (only in this range a venn-table will be calculated)
	if ((K >= 2) & (K <= 4)) {
		#create new lists with the gene-names as list-entries (currently only the distances are in the columns because the summary-table is used)
		newlists <- list()
		for (x in 2:(K + 1)) {
			currentlist <- c()
			for (y in 1:length(genetable[,1])) {
				#check if gene-symbol is in the current list (if the distance is NA, the gene is not in the current list)
				if (!is.na(genetable[,x][y])) {
					currentlist <- append(currentlist, genetable[,1][y], after = y)
				}
			}
			newlists[[colnames(genetable)[x]]] <- currentlist
		}
		
		if(K == 2) {
			#create venn-table for two input lists
			venntable <- matrix(ncol = 2, nrow = 1)
			colnames(venntable) <- c("intersection", "objects")
			venntable[1,1] <- paste(names(newlists)[1], names(newlists)[2], sep = "_")
			venntable[1,2] <- toString(intersect(newlists[[names(newlists)[1]]], newlists[[names(newlists)[2]]])) #L1_L2
			return(list("vennlists" = newlists, "venntable" = venntable))
		}else if(K == 3) {
			#create venn-table for three input lists
			venntable <- matrix(ncol = 2, nrow = 4)
			temp.L1_L2 <- intersect(newlists[[names(newlists)[1]]], newlists[[names(newlists)[2]]]) #to calculate L1_L2_L3
			temp.rownames <- c(paste(names(newlists)[1], names(newlists)[2], names(newlists)[3], sep = "_"),
					paste(names(newlists)[1], names(newlists)[2], sep = "_"),
					paste(names(newlists)[1], names(newlists)[3], sep = "_"),
					paste(names(newlists)[2], names(newlists)[3], sep = "_"))
			venntable[,1] <- temp.rownames
			colnames(venntable) <- c("intersection", "objects")
			venntable[1,2] <- toString(intersect(temp.L1_L2, newlists[[names(newlists)[3]]])) #L1_L2_L3
			venntable[2,2] <- toString(intersect(newlists[[names(newlists)[1]]], newlists[[names(newlists)[2]]])) #L1_L2
			venntable[3,2] <- toString(intersect(newlists[[names(newlists)[1]]], newlists[[names(newlists)[3]]])) #L1_L3
			venntable[4,2] <- toString(intersect(newlists[[names(newlists)[2]]], newlists[[names(newlists)[3]]])) #L2_L3
			return(list("vennlists" = newlists, "venntable" = venntable))
		}else if(K == 4) {
			#create venn-table for four input lists
			venntable <- matrix(ncol = 2, nrow = 9)
			temp.L1_L2 <- intersect(newlists[[names(newlists)[1]]], newlists[[names(newlists)[2]]])
			temp.L3_L4 <- intersect(newlists[[names(newlists)[3]]], newlists[[names(newlists)[4]]])
			temp.rownames <- c(paste(names(newlists)[1], names(newlists)[2], names(newlists)[3], names(newlists)[4], sep = "_"),
					paste(names(newlists)[1], names(newlists)[2], names(newlists)[3], sep = "_"),
					paste(names(newlists)[1], names(newlists)[2], names(newlists)[4], sep = "_"),
					paste(names(newlists)[1], names(newlists)[3], names(newlists)[4], sep = "_"),
					paste(names(newlists)[2], names(newlists)[3], names(newlists)[4], sep = "_"),
					paste(names(newlists)[1], names(newlists)[2], sep = "_"),
					paste(names(newlists)[2], names(newlists)[3], sep = "_"),
					paste(names(newlists)[3], names(newlists)[4], sep = "_"),
					paste(names(newlists)[1], names(newlists)[4], sep = "_"))
			venntable[,1] <- temp.rownames
			colnames(venntable) <- c("intersection", "objects")
			venntable[1,2] <- toString(intersect(temp.L1_L2, temp.L3_L4)) #L1_L2_L3_L4
			venntable[2,2] <- toString(intersect(temp.L1_L2, newlists[[names(newlists)[3]]])) #L1_L2_L3
			venntable[3,2] <- toString(intersect(temp.L1_L2, newlists[[names(newlists)[4]]])) #L1_L2_L4
			venntable[4,2] <- toString(intersect(newlists[[names(newlists)[1]]], temp.L3_L4)) #L1_L3_L4
			venntable[5,2] <- toString(intersect(newlists[[names(newlists)[2]]], temp.L3_L4)) #L2_L3_L4
			venntable[6,2] <- toString(intersect(newlists[[names(newlists)[1]]], newlists[[names(newlists)[2]]])) #L1_L2
			venntable[7,2] <- toString(intersect(newlists[[names(newlists)[2]]], newlists[[names(newlists)[3]]])) #L2_L3
			venntable[8,2] <- toString(intersect(newlists[[names(newlists)[3]]], newlists[[names(newlists)[4]]])) #L3_L4
			venntable[9,2] <- toString(intersect(newlists[[names(newlists)[1]]], newlists[[names(newlists)[4]]])) #L1_L4
			return(list("vennlists" = newlists, "venntable" = venntable))
		}			
	} else {
		return(list("vennlists" = NA, "venntable" = NA))
	}
}





.determineDeltaRange <- function(lists, N, K, v) {
	start.delta <- 0
	stop.delta <- 0
	start.delta.found <- FALSE
	
	#iterate through all possible delta-values
	for (d in 1:N) {
		#check if truncated lists are calculatable with current delta
		#if not calculatable, then the current delta value cannot be applied on the lists (aggmap cannot be drawn)
		if (!.isCalculatable(lists, d, K, v)) {
			if (start.delta.found) {
				return(c(start.delta, d-1))
			}
		} else {
			if (!start.delta.found) {
				start.delta <- d
				start.delta.found <- TRUE
			} else {
				stop.delta <- d
			}
		}
	}
	
	return(c(start.delta, stop.delta))
}





##### added variable lists  - the nontruncated matrix of gene lists, in order to be able to derive the length for color range

.drawAggmap <- function(truncated.lists, N, K, d, v, threshold, lists) {

#delta_symbol = substitute(delta)
#nu_symbol = expression(nu)


if (!is.null(truncated.lists)) #####EBcorr
{
	require(grid)
	wid <- 1
	hei <- 1
	#create list of background coulours to shade the distant-polygons

	

	##### EBcorrection: The color range should represent complete range of the non-truncated list, to give an idea on whether the genes are present on TOP of the whole list, or not. Therefore the range must be the maximal length of non truncated lists:###
	red.white.green <- colorRampPalette(c("red", "orange", "yellow", "white"))(nrow(lists) + 1) #####EBcorr
	
	#get the necessary data
	compared.lists <- truncated.lists$comparedLists
	info <- truncated.lists$info
	grayshaded.lists <- truncated.lists$grayshadedLists
	max.length <- min(length(compared.lists[[1]]), N) #####EBcorr
	
	grid.newpage()
	
	#output the header (the parameters and the colour-gradient)
	pushViewport(viewport(layout = grid.layout(2, 1, widths = c(1, 1), heights = c(0.1, 0.9))))
	pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
	pushViewport(viewport(layout = grid.layout(1, 1, widths = c(1, 1), heights = rep(1, 2))))
	grid.text(substitute(paste(hat(k)[max],"= ", kmax,", N = ", N, ", L = ", K, sep = ""), list(N=N, K=K, kmax=length(compared.lists[[1]]))), x = 0.2, y = 0.65)#####EBcorr
	grid.text(substitute(paste(delta, "= ", a, ", ",nu," = ", b, ", threshold = ", f, "%", sep = ""), list(a=d, b=v, f=threshold)), x = 0.2, y = 0.25)
	upViewport()
	pushViewport(viewport(layout = grid.layout(2, 2, widths = c(1, 1), heights = rep(1, 2))))
	pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
	pushViewport(plotViewport(margins = c(1.5, 1, 1, 1)))
	grid.points(seq(0, 1, length =  nrow(lists)), rep(1, nrow(lists)), pch = 15, gp = gpar(col = red.white.green))#####EBcorr
	pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))	
	grid.text("close", x = 0, y = 0.5, gp = gpar(fontsize = 8, col = "black"))
	grid.text("RANK POSITION", x = 0.5, y = 0.5, gp = gpar(fontsize = 7, col = "black"))
	grid.text("distant", x = 1, y = 0.5, gp = gpar(fontsize = 8, col = "black"))
	upViewport()
	popViewport()
	popViewport()
	upViewport()
	upViewport()
	
	#set the layout of the 
	pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))	
	pushViewport(viewport(layout = grid.layout(max.length + 1, 1+truncated.lists$Ntoplot, widths = rep(1, 1+truncated.lists$Ntoplot) * (max.length), heights = rep(1, 1 + truncated.lists$Ntoplot) * (max.length))))#####EBcorr
	
	#output row-numbers
	for (g in c(1:max.length)) {
		pushViewport(viewport(width = wid, height = hei, layout.pos.col = 1, layout.pos.row = g + 1))
		grid.text(g, gp = gpar(fontsize = 8, col = "black"))
		popViewport()
	}
	
	h.sizes <- c(K:2)
	
	#output the lists
	for (i in 1:length(compared.lists)) {
		#get the name of the list and output it
		listname <- info[2, i]
		pushViewport(viewport(width = wid, height = hei, layout.pos.col = i + 1, layout.pos.row = 1))
		grid.text(listname, gp = gpar(fontsize = 8, col = "black"))
		popViewport()
		
		#output the list (check if to output a reference or truncated list)
		if (info[3,i] == "R") {
			#output reference list
			for (y in 1:min(length(compared.lists[[i]]), max.length)) {
				pushViewport(viewport(width = wid, height = hei, layout.pos.col = i + 1, layout.pos.row = y + 1))				
				#check if the gene-symbol has tho be gray-shaded
				bg <- "white"
				if (!is.na(grayshaded.lists[[i]][y])) { if(grayshaded.lists[[i]][y]==1) {
					bg <- "grey"
				}}
				grid.rect(gp = gpar(col = "black", fill = bg))
				grid.text(compared.lists[[i]][y], gp = gpar(fontsize = 8, col = "black", fontface = "bold"))
				popViewport()
			}
		} else if (info[3,i] == "T") {
			#output truncated list
			for (y in 1:min(length(compared.lists[[i]]), max.length)) {
				distance <- compared.lists[[i]][y]
				pushViewport(viewport(width = wid, height = hei, layout.pos.col = i + 1, layout.pos.row = y + 1))
				if (!is.na(distance)) {
					grid.polygon(x = c(1, 1, 0), y = c(0, 1, 0), gp = gpar(col = "black", fill = red.white.green[distance + 1]))
					grid.text(distance, x = 0.9, y = 0.4, gp = gpar(cex = 0.7))
				} else {
					grid.polygon(x = c(1, 1, 0), y = c(0, 1, 0), gp = gpar(col = "black", fill = "white"))
					grid.text("NA", x = 0.9, y = 0.4, gp = gpar(cex = 0.7))
				}
				bg <- "white"
				if (!is.na(grayshaded.lists[[i]][y])) { if(grayshaded.lists[[i]][y]==1) {
					bg <- "grey"
				}}
				grid.polygon(x = c(0, 1,0), y = c(1, 1, 0), gp = gpar(col = "black", fill = bg))
				popViewport()				
			}
		}
	}
	upViewport()
	upViewport()
   } else {grid.newpage()
	red.white.green <- colorRampPalette(c("red", "orange", "yellow", "white"))(nrow(lists) + 1) #####EBcorr
	#output the header (the parameters and the colour-gradient)
	pushViewport(viewport(layout = grid.layout(2, 1, widths = c(1, 1), heights = c(0.1, 0.9))))
	pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
	pushViewport(viewport(layout = grid.layout(1, 1, widths = c(1, 1), heights = rep(1, 2))))
	grid.text(substitute(paste(k[max],"= NA, N = ", N, ", L = ", K, sep = ""), list(N=N, K=K)), x = 0.2, y = 0.65)#####EBcorr
	grid.text(substitute(paste(delta, "= ", a, ", ",nu," = ", b, ", threshold = ", f, "%", sep = ""), list(a=d, b=v, f=threshold)), x = 0.2, y = 0.3)
	upViewport()
	pushViewport(viewport(layout = grid.layout(2, 2, widths = c(1, 1), heights = rep(1, 2))))
	pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
	pushViewport(plotViewport(margins = c(1.5, 1, 1, 1)))
	grid.points(seq(0, 1, length =  nrow(lists)), rep(1, nrow(lists)), pch = 15, gp = gpar(col = red.white.green))#####EBcorr
	pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))	
	grid.text("close", x = 0, y = 0.5, gp = gpar(fontsize = 8, col = "black"))
	grid.text("RANK POSITION", x = 0.5, y = 0.5, gp = gpar(fontsize = 7, col = "black"))
	grid.text("distant", x = 1, y = 0.5, gp = gpar(fontsize = 8, col = "black"))
	upViewport()
	popViewport()
	popViewport()
	upViewport()
	upViewport()
	#set the layout of the 
	pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))	
	pushViewport(viewport(layout = grid.layout(1, 2, widths = c(1), heights = c(1,1))))#####EBcorr
	grid.text(substitute(paste("On selected ", delta ,", no overlap found", sep="")))
	upViewport()
	upViewport()
	upViewport()  } #end for is.null #####EBcorr
}


.isCalculatable <- function(lists, d, K, v) {
	res.j0.temp <- j0.multi(lists, d, v)
	res.temp <- as.matrix(res.j0.temp$K)
	
	temp2 = c()
	#calculate the reference-lists
	for (i in c(1:K)) {
		if (i > 1) {
			res.temp.temp <- res.temp[-c(which(!is.na(match(res.temp[,1], names(temp2)))), which(!is.na(match(res.temp[,2], names(temp2))))),]
		} else {
			res.temp.temp <- res.temp
		}
		temp <- tapply(as.numeric(res.temp.temp[,4]), res.temp.temp[,1], function(x) max(x, na.rm = TRUE))
		which.max.temp <- which.max(temp)
		temp2 <- c(temp2, temp[which.max.temp])
	}
	
	lists.to.remove <- ""
	truncation.points <- c()
	
	#calculate the truncated-lists for each block (for each reference-list)
	for (first.list.name in names(temp2)) {
		res.temp.selection <- res.temp[which(res.temp[,1] == first.list.name), ]
		res.temp.selection <- res.temp.selection[which(is.na(match(res.temp.selection[,2], lists.to.remove))),]
		if (is.matrix(res.temp.selection)) {
			ranked.list.names <- c(first.list.name, res.temp.selection[order(res.temp.selection[,4], decreasing = TRUE),2])
			lists.rank <- rank(ranked.list.names)
		} else {
			ranked.list.names <- c(first.list.name, res.temp.selection[2])
			lists.rank <- rank(ranked.list.names)
		}
		
		#get the truncation-points for the truncated lists in the current block
		for(currennt.truncated.list in ranked.list.names[2:length(ranked.list.names)]) {
			truncation.points <- c(truncation.points, as.numeric(res.temp[res.temp[,1] == first.list.name & res.temp[,2] == currennt.truncated.list,4]))
		}
		
		lists.to.remove <- c(lists.to.remove, first.list.name)
	}
	
	#check if a truncation-point is missing - if true, the delta value cannot be applied on the lists (aggmap cannot be drawn)
	if(NA %in% truncation.points) {
		return(FALSE)
	}else {
		return(TRUE)
	}
}






## OLD DATA


## aggmap<-function(lists, res.j0.temp, N, K, delta, v)
## {

## wid=1
## hei=1
## res.temp =  as.matrix(res.j0.temp$K)
## max.j0.est = res.j0.temp$maxK

## local.max.j0.est = rep(K, K)
## names(local.max.j0.est) = names(lists)
## temp2=c()

## for (i in c(1:K))
## {
## if(i>1){res.temp.temp = res.temp[-c(which(!is.na(match(res.temp[,1],names(temp2)))),which(!is.na(match(res.temp[,2],names(temp2))))),]}else{res.temp.temp=res.temp}
## temp = tapply(as.numeric(res.temp.temp[,4]), res.temp.temp[,1], function(x) max(x, na.rm=TRUE))
## which.max.temp = which.max(temp)
## temp2 = c(temp2, temp[which.max.temp])
## }

## local.max.j0.est[names(temp2)]=order(temp2, decreasing=T)
## first.order.local.max.j0.est = order(local.max.j0.est)

## red.white.green = colorRampPalette(c("red", "orange","yellow","white"))(N+1)
## library(grid)
## grid.newpage()

## # defining basic layout with two rows and one column
## pushViewport(viewport(layout=grid.layout(2, 1, widths=c(1,1),heights=c(0.1,0.9))))

## ##first row - defining header (delta, nu and color key level) 
## 	pushViewport(viewport(layout.pos.row=1, layout.pos.col=1))
## ### first column - delta and nu 
## 		pushViewport(viewport(layout=grid.layout(1, 1, widths=c(1,1),heights=rep(1,2))))
## 		grid.text(paste("delta = ",delta,",  v =",v,sep=""), x=0.2, y=0.5)
## 		upViewport()		

## ### second column - color key
## 		pushViewport(viewport(layout=grid.layout(1, 2, widths=c(1,1),heights=rep(1,2)))) # defines two columns in the second column
## 			pushViewport(viewport(layout.pos.row=1, layout.pos.col=2)) # fills the second column with the key
## 				pushViewport(plotViewport(margins=c(1.5,1,1,1)))
## 				grid.xaxis(at=seq(0, 1, length=11), gp=gpar(cex=0.5,tck=0.1, line=0.1))
## 				grid.points(seq(0, 1, length=N), rep(1, N), pch=15, gp=gpar(col=red.white.green))
## 				popViewport()
## 			popViewport()
## 		upViewport()
## 	upViewport() #coming back to the basic level

## ## second row - plotting genes
## 	pushViewport(viewport(layout.pos.row=2, layout.pos.col=1))

## ### defining the layout in second row
## 		pushViewport(viewport(layout=grid.layout(max.j0.est+v+1, 1+factorial(K)/2+2, widths=rep(1,1+factorial(K)/2+2)*(max.j0.est+v+1),heights=rep(1,1+factorial(K)/2+2)*(max.j0.est+v+1))))

## 			for (g in c(1:(max.j0.est+v)))
## 				{
## 				pushViewport(viewport(width=wid, height=hei, layout.pos.col=1, layout.pos.row=g+1))
## 				grid.text(g, gp=gpar(fontsize=8, col="black"))
## 				popViewport()
## 				}# end for g
## 		h.sizes=c(K:2)
## 		l=1
## 		lists.to.remove = ""

## 	for (first.list.name in names(temp2))
## 		{
## 		# selecting the order of pairs of lists to the first one
## 		res.temp.selection = res.temp[which(res.temp[,1]==first.list.name),]
## 		res.temp.selection = res.temp.selection[which(is.na(match(res.temp.selection[,2],lists.to.remove))),]
## 			if(is.matrix(res.temp.selection)){
## 			ranked.list.names = c(first.list.name, res.temp.selection[order(res.temp.selection[,4], decreasing=TRUE),2])
## 			lists.rank = rank(ranked.list.names)
## 			order.local.max.j0.est = order(as.numeric(res.temp.selection[,4]), decreasing=TRUE)
## 			}else{
## 			ranked.list.names = c(first.list.name,res.temp.selection[2])
## 			lists.rank = rank(ranked.list.names)
## 			order.local.max.j0.est = 1
## 			}

## 		gene.names.to.plot = lists[,first.list.name][1:(temp2[first.list.name]+v)]


## 		### plotting gene names ###
## 		pushViewport(viewport(width=wid, height=hei, layout.pos.col=l+1,layout.pos.row=1))
## 		grid.text(ranked.list.names[1], gp=gpar(fontsize=8, col="black"))
## 		popViewport()

## 		for (g in c(1:(length(gene.names.to.plot))))
## 				{
## 				pushViewport(viewport(width=wid, height=hei, layout.pos.col=l+1, layout.pos.row=g+1))
## 				grid.rect(gp=gpar(col="black"))
## 				grid.text(gene.names.to.plot[g], gp=gpar(fontsize=8, col="black"))
## 				popViewport()
## 				}# end for g
		
## 		for (k in c(2:K))
## 		{
## 		pushViewport(viewport(width=wid, height=hei, layout.pos.col=l+k, layout.pos.row=1))
## 		grid.text(ranked.list.names[k], gp=gpar(fontsize=8, col="black"))
## 		popViewport()
## 		n.genes.to.plot = as.numeric(res.temp[res.temp[,1]==first.list.name & res.temp[,2]==ranked.list.names[k],4])+v
## 			for (g in c(1:n.genes.to.plot))
## 			{
## 			j = match(gene.names.to.plot[g],lists[,ranked.list.names[k]])
## 			distance = abs(g-j)
## 				if (!is.na(j) & j<(n.genes.to.plot+v)) # if gene is in the second list AND is in the top j0+v genes of the second list
## 					{i="grey"}else{i="white"}
				
## 				pushViewport(viewport(width=wid, height=hei, layout.pos.col=l+k, layout.pos.row=g+1))
## 				if (length(distance)>0)
## 					{
## 					grid.polygon(x=c(0,1,1),y=c(1,1,0), gp=gpar(col="black", fill=red.white.green[distance+1]))#upper right triangle
## 					grid.text(distance,x=0.8, y=0.6, gp=gpar(cex=0.7))
## 					}else
## 					{
## 					grid.polygon(x=c(0,1,1),y=c(1,1,0), gp=gpar(col="black", fill="white"))#upper right triangle
## 					grid.text("NA",x=0.8, y=0.6, gp=gpar(cex=0.7))
## 					}
## 				grid.polygon(x=c(0,0,1),y=c(1,0,0), gp=gpar(col="black", fill=i))#bottom left triangle
## 				popViewport()
## 			}# end for g
## 		}# end for k
	
## 			l=l+K
## 			K=K-1
## 			lists.to.remove = c(lists.to.remove, first.list.name)

## 		}# end for first list name


## 			upViewport()
		
## 	upViewport()
## } # end of the function
