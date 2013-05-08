TopKListsGUI <- function(lists, autorange.delta = FALSE, override.errors = FALSE, aggmap.pdf.size = c(9, 8), venndiag.pdf.size = c(7, 7), venndiag.size = c(380, 420), aggmap.size = c(870, 440), gui.size = c(900, 810), directory = NULL, venndiag.res = 70, aggmap.res = 100, maxd = 500) {
  options("guiToolkit"="RGtk2")

  ##setting up the directory
  if(is.null(directory)) {
    directory <- paste(getwd(), "/TopKLists-temp",sep="")
    if(!file.exists(directory)) dir.create(directory)
  }
  message(paste("Writing files to", directory))

  
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
  current.arguments <- environment() #contains all the arguments (N,L,v,....) of the current calculation

                                        #create input area for arguments
  calcframe <- gframe("ARGUMENTS", pos = 0.5, container = maingroup)
  calcgroup <- ggroup(horizontal = FALSE, container = calcframe, expand = TRUE)
  agroup <- ggroup(horizontal = TRUE, container = calcgroup, expand = TRUE)
  a1frame <- gframe("Calculated from given list", container = agroup, expand = TRUE)
  a2frame <- gframe("Choose range for delta (depends on nu)", container = agroup, expand = TRUE)
  a3frame <- gframe("Threshold", container = agroup, expand = TRUE)

  wid$calculate <- gbutton("Calculate", container = calcgroup, handler = function(h,...) { calculate.all.truncationlists() })	
  svalue(calcgroup) <- 10

                                        #layout of the first argument input area
  tbl1 <- glayout(spacing = 10,container = a1frame)
  tbl1[2,2, anchor = c(1,0)] <- glabel("N", container = tbl1)
  tbl1[2,3] <- wid$N <- gspinbutton(from = 0, to = dim(lists)[1], by = 1, value = dim(lists)[1], container = tbl1)
  tbl1[2,4, anchor = c(-1,0)] <- glabel("maximal number of objects" ,container = tbl1)
  tbl1[3,2, anchor = c(1,0)] <- glabel("L", container = tbl1)
  tbl1[3,3] <- wid$L <- gspinbutton(from = 0, to = dim(lists)[2], by = 1, value = dim(lists)[2], container = tbl1)
  tbl1[3,4, anchor = c(-1,0)] <- glabel("number of lists", container = tbl1)

                                        #layout of the second argument input area
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

                                        #layout of the third argument input area
  tbl3 <- glayout(spacing = 10, container = a3frame)
  tbl3[2,2, anchor = c(1,0)] <- glabel("min", container = tbl3)
  tbl3[2,3] <- wid$threshold <- gedit(text = "50", width = 3, container = tbl3)
  tbl3[2,4, anchor = c(-1,0)] <- glabel("%" ,container = tbl3)
  tbl3[3,2:4, anchor = c(1,0)] <- glabel("to gray-shade gene symbol", container = tbl3)

                                        #create the delta slider for choosing the desired delta value (and to show the calculated data for the selected delta value)
  sldrframe <- gframe("DELTA-SLIDER", pos = 0.5, container = maingroup)
                                        #the parameters 'from' and 'to' are dummy values because the slider will be updated after the calculation
  wid$delta.slider <- gslider(from = 1, to = 2, by = 1, container = sldrframe, expand = TRUE, handler = function(h,...) { load.data() })
  enabled(wid$delta.slider) <- FALSE

                                        #create four tabs for presenting the results of the calculation
  nb <- gnotebook(container = maingroup, expand = TRUE)
  aggmapg <- ggroup(horizontal = FALSE, container = nb, label = "Aggregation map")
  summtblg <- ggroup(horizontal = FALSE, container = nb, label = "Summary table")
  venng <- ggroup(horizontal = TRUE, container = nb, label = "Venn-diagram & Venn-table")
  svalue(nb) <- 1 #set the first tab as the selected tab

                                        #create a status bar to inform the user of the current status
  wid$status <- gtkStatusbar()
  add(maingroup, wid$status)
  info <- wid$status$getContextId("info")
  wid$status$push(info, "GUI started")

                                        #create a progress bar to show the progress of the calculation
  wid$progress <- gtkProgressBar(show = TRUE)
  add(maingroup, wid$progress)	

                                        #update the allowed range for delta (based on the v-value)
  update.deltarange <- function() {
                                        #check if allowed range for delta should be updated
    if (autorange.delta) {
      deltarange <- .determineDeltaRange(lists, as.numeric(svalue(wid$N)), as.numeric(svalue(wid$L)), as.numeric(svalue(wid$v)))
      svalue(wid$delta.start) <- deltarange[1]
      svalue(wid$delta.stop) <- deltarange[2]
      enabled(wid$delta.start) <- TRUE
      enabled(wid$delta.stop) <- TRUE
      wid$status$push(info, expression(paste("Recalculated allowed range for delta with ",nu," =", ss), list(ss=svalue(wid$v))))
    }
  }

                                        #update the calculated data set (aggmap, summary-table and Venn) in the three tabs
  load.data <- function(truncated.lists) {		
                                        #update is only performed when delta slider is enabled
    if (enabled(wid$delta.slider)) {
                                        #check if calculated data set can be shown in the GUI (no file = error in the calculation for this delta-value)
      if (file.exists(paste(directory, "/N" , current.arguments$val.N, "_L", current.arguments$val.L, "_delta", svalue(wid$delta.slider), "_v", current.arguments$val.v, "_thrshld", current.arguments$val.threshold, ".Rdata", sep = ""))) {
                                        #load the calculated data set
        load(paste(directory, "/N", current.arguments$val.N, "_L", current.arguments$val.L, "_delta", svalue(wid$delta.slider), "_v", current.arguments$val.v, "_thrshld", current.arguments$val.threshold, ".Rdata", sep = ""))

                                        #first tab: show aggmap
        if (!is.null(wid$aggmap)) {
					#delete current aggmap image
          delete(aggmapg, wid$aggmap)
        }
        wid$aggmap <- gimage(paste(directory, "/aggmap_N", current.arguments$val.N, "_L", current.arguments$val.L, "_delta", svalue(wid$delta.slider), "_v", current.arguments$val.v, "_thrshld", current.arguments$val.threshold, ".png", sep = ""), container = aggmapg)


                                       #second tab: show summary-table
        if (!is.null(wid$sumtable)) {
					#delete current summary-table
          delete(summtblg, wid$sumtable)
        }
        wid$sumtable <- gtable(truncated.lists$summarytable, container = summtblg, expand = TRUE)

                                        #third tab: show Venn-diagram and Venn-table (only if L <= 4)
        if (current.arguments$val.L <= 4) {			
          if (!is.null(wid$venndiag) & !is.null(wid$venntbl)) {
                                        #delete current Venn-diagram and Venn-table
            delete(venng, wid$venndiag)
            delete(venng, wid$venntbl)
          }
          if (!is.null(wid$nodraw) & !is.null(wid$nodrawimage)) {
                                        #delete the icon and message
            delete(venng, wid$nodrawimage)
            delete(venng, wid$nodraw)
          }
          wid$venndiag <- gimage(paste(directory, "/venn_N", current.arguments$val.N, "_L", current.arguments$val.L, "_delta", svalue(wid$delta.slider), "_v", current.arguments$val.v, "_thrshld", current.arguments$val.threshold, ".png", sep = ""), container = venng)
          wid$venntbl <- gtable(truncated.lists$venntable, container = venng, expand = TRUE)

        } else {
					#output message that Venn-diagram and Venn-table are not shown for L > 4
          if (is.null(wid$nodrawimage) & is.null(wid$nodraw)) {
            wid$nodrawimage <- gimage("info", dirname = "stock", container = venng)
            wid$nodraw <- glabel("There are more than four input lists.\nVenn-diagram and Venn-table are only drawn to a maximum of four input lists.", container = venng)
          }
        }# end for third tab if
  
      } else {
        gmessage("There is no data to show in the GUI!\nThis is caused by an error in the calculation.\nYou may have to adjust the arguments!", title = "Loading data failed", icon = "error")
      }
    }
  }
  
                                        #calculates aggmap, summary-table and Venn for each delta in the given range for delta
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

                                        #calculate data set for each delta in the given range for delta
    temp.tocalc <- (as.numeric(svalue(wid$delta.stop)) - as.numeric(svalue(wid$delta.start))) + 1
    temp.loopcounter <- 0
    wid$error <- FALSE
    for (i in as.numeric(svalue(wid$delta.start)):as.numeric(svalue(wid$delta.stop))) {
      gtkMainIterationDo(FALSE)
                                        #try to calculate data set for current delta
      tryCatch({
        truncated.lists <- calculateDataSet(lists, as.numeric(svalue(wid$L)), i, as.numeric(svalue(wid$v)), as.numeric(svalue(wid$threshold)))
                                        #save calculated truncated lists for current delta in the destination directory
        save(truncated.lists, file = paste(directory, "/N", as.numeric(svalue(wid$N)), "_L", as.numeric(svalue(wid$L)), "_delta", i, "_v", as.numeric(svalue(wid$v)), "_thrshld", as.numeric(svalue(wid$threshold)), ".Rdata", sep = ""))
                                        #draw the aggmap and save it as pdf-image in the destination directory
        pdf(file = paste(directory, "/aggmap_N", as.numeric(svalue(wid$N)), "_L", as.numeric(svalue(wid$L)), "_delta", i, "_v", as.numeric(svalue(wid$v)), "_thrshld", as.numeric(svalue(wid$threshold)), ".pdf", sep = ""), width = aggmap.pdf.size[1], height = aggmap.pdf.size[2])
        aggmap(truncated.lists, N=as.numeric(svalue(wid$N)), L=as.numeric(svalue(wid$L)), d=i, v=as.numeric(svalue(wid$v)), threshold=as.numeric(svalue(wid$threshold)), lists)
        dev.off()

 	  png(filename = paste(directory, "/aggmap_N", as.numeric(svalue(wid$N)), "_L", as.numeric(svalue(wid$L)), "_delta", i, "_v", as.numeric(svalue(wid$v)), "_thrshld", as.numeric(svalue(wid$threshold)), ".png", sep = ""), width = aggmap.size[1], height = aggmap.size[2], res=aggmap.res)
        aggmap(truncated.lists, N=as.numeric(svalue(wid$N)), L=as.numeric(svalue(wid$L)), d=i, v=as.numeric(svalue(wid$v)), threshold=as.numeric(svalue(wid$threshold)), lists)
        dev.off()

                                      #draw the Delta-plot and save it as pdf-image in the destination directory
        pdf(file = paste(directory, "/deltaplot_N", as.numeric(svalue(wid$N)), "_L", as.numeric(svalue(wid$L)), ".pdf", sep = ""), width = aggmap.pdf.size[1], height = aggmap.pdf.size[2])
        Mdelta=deltaplot(lists)
        dev.off()

        save(Mdelta, file=paste(directory, "/Mdelta.rdata", sep=""))
        for (aname in names(Mdelta))
        {
	write.table(Mdelta[[aname]], file=paste(directory, "/Mdelta",aname,".txt", sep=""), sep="\t", row.names=FALSE)
        }

 	png(filename = paste(directory, "/deltaplot_N", as.numeric(svalue(wid$N)), "_L", as.numeric(svalue(wid$L)), ".png", sep = ""), width = aggmap.size[1], height = aggmap.size[2], res=aggmap.res)
        Mdelta=deltaplot(lists)
        dev.off()
        rm(Mdelta)

                                        #draw the Venn-diagram and save it as pdf-image in the destination directory (only if L is between 2 and 4 & if there is a j0 estimated)
        if((as.numeric(svalue(wid$L)) >= 2) & (as.numeric(svalue(wid$L)) <= 4) & !is.null(truncated.lists)) {
          pdf(file = paste(directory, "/venn_N", as.numeric(svalue(wid$N)), "_L", as.numeric(svalue(wid$L)), "_delta", i, "_v", as.numeric(svalue(wid$v)), "_thrshld", as.numeric(svalue(wid$threshold)), ".pdf", sep = ""), width = venndiag.pdf.size[1], height = venndiag.pdf.size[2], bg = "transparent")
          venn(truncated.lists$vennlists)
          dev.off()
          png(filename = paste(directory, "/venn_N", as.numeric(svalue(wid$N)), "_L", as.numeric(svalue(wid$L)), "_delta", i, "_v", as.numeric(svalue(wid$v)), "_thrshld", as.numeric(svalue(wid$threshold)), ".png", sep = ""), width = venndiag.size[1], height = venndiag.size[2], bg = "transparent", res=venndiag.res)
          venn(truncated.lists$vennlists)
          dev.off()
        }else{
          pdf(file = paste(directory, "/venn_N", as.numeric(svalue(wid$N)), "_L", as.numeric(svalue(wid$L)), "_delta", i, "_v", as.numeric(svalue(wid$v)), "_thrshld", as.numeric(svalue(wid$threshold)), ".pdf", sep = ""), width = venndiag.pdf.size[1], height = venndiag.pdf.size[2], bg = "transparent")
          plot(1,1, type="n", axes=FALSE)
          text(1,1,"No overlap on selected parameters")
          dev.off()

 	 png(filename = paste(directory, "/venn_N", as.numeric(svalue(wid$N)), "_L", as.numeric(svalue(wid$L)), "_delta", i, "_v", as.numeric(svalue(wid$v)), "_thrshld", as.numeric(svalue(wid$threshold)), ".png", sep = ""), width = venndiag.size[1], height = venndiag.size[2], bg = "transparent", res=venndiag.res)
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

                                        #update progress bar
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
                                        #(this prevents errors when sliding through the delta-values - otherwise changing the arguments in the GUI could cause errors when showing the calculated data)
    current.arguments$val.N <- as.numeric(svalue(wid$N))
    current.arguments$val.L <- as.numeric(svalue(wid$L))
    current.arguments$val.v <- as.numeric(svalue(wid$v))
    current.arguments$val.threshold <- as.numeric(svalue(wid$threshold))
    
                                        #show the data in the tabs
    load.data()
  }#end calculate.all.truncationlists
  
                                        #show the GUI
  visible(win) <- TRUE
}




.calculateVennValues <- function(genetable, L) {
                                        #check if L is between 2 and 4 (only in this range a Venn-table will be calculated)
  if ((L >= 2) & (L <= 4)) {
                                        #create new lists with the object-symbols or names as list-entries (currently only the distances are in the columns because the summary-table is used)
    newlists <- list()
    for (x in 2:(L + 1)) {
      currentlist <- c()
      for (y in 1:length(genetable[,1])) {
                                        #check if object-symbol is in the current list (if the distance is NA, the object is not in the current list)
        if (!is.na(genetable[,x][y])) {
          currentlist <- append(currentlist, genetable[,1][y], after = y)
        }
      }
      newlists[[colnames(genetable)[x]]] <- currentlist
    }
    
    if(L == 2) {
                                        #create Venn-table for two input lists
      venntable <- matrix(ncol = 2, nrow = 1)
      colnames(venntable) <- c("intersection", "objects")
      venntable[1,1] <- paste(names(newlists)[1], names(newlists)[2], sep = "_")
      venntable[1,2] <- toString(intersect(newlists[[names(newlists)[1]]], newlists[[names(newlists)[2]]])) #L1_L2
      return(list("vennlists" = newlists, "venntable" = venntable))
    }else if(L == 3) {
                                        #create Venn-table for three input lists
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
    }else if(L == 4) {
                                        #create Venn-table for four input lists
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





.determineDeltaRange <- function(lists, N, L, v) {
  start.delta <- 0
  stop.delta <- 0
  start.delta.found <- FALSE
  
                                        #iterate through all possible delta-values
  for (d in 1:N) {
                                        #check if truncated lists are calculatable with current delta
                                        #if not calculatable, then the current delta value cannot be applied on the lists (aggmap cannot be drawn)
    if (!.isCalculatable(lists, d, L, v)) {
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

###function that generates Delta-plot and Delta-matrix
##if subset.plotted is NA no subplots are created
deltaplot<-function(lists, mind=0, maxd=NULL, perc.subplot=50, subset.plotted=NA)
{
  if (is.null(subset.plotted) & nrow(lists) > 200){
    stop("Subset for calculating zero count (subset.plotted) has to be specified \n")
  }
  if (is.null(subset.plotted) & nrow(lists) < 201){
    warning(paste("Subset of the lists for calculating zero count is not specified, using", nrow(lists)), "\n")
    subset.plotted <- nrow(lists)
  }
  if (is.null(maxd)) {
    cat(paste("The maximum for delta not specified, using",nrow(lists)*0.25), "\n")
    maxd=c(nrow(lists)*0.25)
  }
  if (maxd>nrow(lists)) {
    cat(paste("The maximum for delta you specified is larger than the number of objects in your lists. Maxd changed to",nrow(lists)*0.25), "\n")
    maxd=c(nrow(lists)*0.25)
  }
  if(is.null(names(lists)) | any(names(lists)=="")){
    names(lists) <- paste("L",1:ncol(lists),sep="")
    warning(paste("List names not be given or incorrect. Replaced by 'L1', 'L2',... L",ncol(lists),'\n'))
  }

  if(is.null(subset.plotted) | !is.na(subset.plotted)){
    lists = lists[1:subset.plotted,] ## takes only specified subset of input list  
  }
  
  ## calculate zero count for each pair of lists  
  deltas = c(mind:maxd)
  Mdelta = list()
  xxs = list()
  n=ncol(lists)
  a = n*(n-1)
  for (i in 1:(ncol(lists)-1))
    {
      for (j in (i+1):ncol(lists))
        {
			x11()
			par(mfrow=c(1,2))
        ## zero count calculation and deltaplot for LiLj and LjLi (in one window)
			## LiLj
              Mdelta.temp = data.frame(Object=c(as.character(lists[,i]), "#zeros"), L1=c(c(1:nrow(lists)), NA), L2 = c(match(lists[,i],lists[,j]), NA))
              names(Mdelta.temp)[2:3] = c(paste("L",i, sep=""),paste("L",j, sep=""))
              xx = c()
              for (d in deltas)
                {	
                  a = prepareIdata(lists[,c(i,j)],d=d)
                  x = table(as.numeric(a$Idata))['0']
                  xx = c(xx,x)
                  Mdelta.temp[,paste("delta_",d)] = c(a$Idata, x)
                }# end for d
			  Mdelta[[paste("L",i,"L",j, sep="")]] = Mdelta.temp
              xxs[[paste("L",i,"L",j, sep="")]] = xx  ##saving xx for plotting single deltaplot with subplot in the corner
	    	  par(mar=c(5,5,1,1))
              plot(deltas,xx, xlab=expression(delta), ylab="# of 0's", las=1,cex.axis=0.7, main=paste("L",i, " vs L",j, sep=""))

			## LjLi  
              Mdelta.temp = data.frame(Object=c(as.character(lists[,j]), "#zeros"), L1=c(c(1:nrow(lists)), NA), L2 = c(match(lists[,j],lists[,i]), NA))
              names(Mdelta.temp)[2:3] = c(paste("L",j, sep=""),paste("L",i, sep=""))
              xx = c()
              for (d in deltas)
                {	
                  a = prepareIdata(lists[,c(j,i)],d=d)
                  x = table(a$Idata)['0']
                  xx = c(xx,x)
                  Mdelta.temp[,paste("delta_",d)] = c(a$Idata, x)
                }# end for d
			  Mdelta[[paste("L",j,"L",i, sep="")]] = Mdelta.temp
              xxs[[paste("L",j,"L",i, sep="")]] = xx  ##saving xx for plotting single deltaplot with subplot in the corner
	    	  par(mar=c(5,5,1,1))
              plot(deltas,xx, xlab=expression(delta), ylab="# of 0's", las=1,cex.axis=0.7, main=paste("L",j, " vs L",i, sep=""))		  
        }# end for j
    }# end for i

	## plot deltaplot
	par(mfrow=c(1,2))
    	
	
  if(!is.na(subset.plotted)){
    ## deltaplot with subplot in the top right corner:
    for (i in 1:ncol(lists))
      {
        for (j in 1:ncol(lists))
          {
            if (i!=j){
              x11()
              par(mar=c(5,5,1,1))
              plot(deltas,xxs[[paste("L",i,"L",j,sep="")]], xlab=expression(delta), ylab="# of 0's", las=1,cex.axis=0.7, main=paste("L",i, " vs L",j, sep=""))			
              extremes = par("usr")
              dimen = par("pin")					   
              subplot(plot(deltas[1:((perc.subplot/100)*length(deltas))],xxs[[paste("L",i,"L",j,sep="")]][1:((perc.subplot/100)*length(deltas))], xlab="", ylab="", las=1, cex.axis=0.7) , extremes[2], extremes[4], size = c(dimen[1]*0.5, dimen[2]*0.4),hadj=1, vadj=1, pars=list(col="black", mar=c(5,5,1,1)))   
            }
          }
      }
  }
  return(Mdelta)

}#end deltaplot

## Aggregation map (aggmap) for the graphical representation of more than two top-k lists obtained from the method by Hall and Schimek (2012) - it applies Paul Murrell's grid package
## Input variables 
# L - number of lists
# N - maximal number of objects
# v - value of tuning parameter nu used for j0 estimation
# res.temp - data frame, result from the Hall and Schimek (2012) algorithm for all pairings of lists, 4 columns, the first and second column contain the list names that were compared, the third column contains the selected nu, and the fourth column the estimated j0
# The lists data frame contains the complete set of L ordered input lists - column names are obligatory

aggmap <- function(truncated.lists, N, L, d, v, threshold, lists) {

                                        #delta_symbol = substitute(delta)
                                        #nu_symbol = expression(nu)

  if (!is.null(truncated.lists)) 
    {
      require(grid)
      wid <- 1
      hei <- 1
                                        #create list of background colours to shade the distance-polygons  

      red.white.green <- colorRampPalette(c("red", "orange", "yellow", "white"))(nrow(lists) + 1) 
      
                                        #get the necessary data
      compared.lists <- truncated.lists$comparedLists
      info <- truncated.lists$info
      grayshaded.lists <- truncated.lists$grayshadedLists
      max.length <- min(length(compared.lists[[1]]), N) 

   #  
      grid.newpage()
      
                                        #output the header (the parameters and the colour gradient)
      pushViewport(viewport(layout = grid.layout(2, 1, widths = c(1, 1), heights = c(0.1, 0.9))))
      pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
      pushViewport(viewport(layout = grid.layout(1, 1, widths = c(1, 1), heights = rep(1, 2))))
      grid.text(substitute(paste(hat(k)[max],"= ", kmax,", N = ", N, ", L = ", L, sep = ""), list(N=N, L=L, kmax=length(compared.lists[[1]]))), x = 0.2, y = 0.65)
      grid.text(substitute(paste(delta, "= ", a, ", ",nu," = ", b, ", threshold = ", f, "%", sep = ""), list(a=d, b=v, f=threshold)), x = 0.2, y = 0.25)
      upViewport()
      pushViewport(viewport(layout = grid.layout(2, 2, widths = c(1, 1), heights = rep(1, 2))))
      pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
      pushViewport(plotViewport(margins = c(1.5, 1, 1, 1)))
      grid.points(seq(0, 1, length =  nrow(lists)), rep(1, nrow(lists)), pch = 15, gp = gpar(col = red.white.green))
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
      pushViewport(viewport(layout = grid.layout(max.length + 1, 1+truncated.lists$Ntoplot, widths = rep(1, 1+truncated.lists$Ntoplot) * (max.length), heights = rep(1, 1 + truncated.lists$Ntoplot) * (max.length))))
      
                                        #output row-numbers
      for (g in c(1:max.length)) {
        pushViewport(viewport(width = wid, height = hei, layout.pos.col = 1, layout.pos.row = g + 1))
        grid.text(g, gp = gpar(fontsize = 8, col = "black"))
        popViewport()
      }
      
      h.sizes <- c(L:2)
      
                                        #output the lists
      for (i in 1:length(compared.lists)) {
                                        #get the name of the list and output the list name
        listname <- info[2, i]
        pushViewport(viewport(width = wid, height = hei, layout.pos.col = i + 1, layout.pos.row = 1))
        grid.text(listname, gp = gpar(fontsize = 8, col = "black"))
        popViewport()
        
                                        #output the list (check output for reference or truncated list)
        if (info[3,i] == "R") {
                                        #output reference list
          for (y in 1:min(length(compared.lists[[i]]), max.length)) {
            pushViewport(viewport(width = wid, height = hei, layout.pos.col = i + 1, layout.pos.row = y + 1))				
                                        #check if the object-symbol or name has to be gray-shaded
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
            red.white.green <- colorRampPalette(c("red", "orange", "yellow", "white"))(nrow(lists) + 1) 
                                        #output the header (the parameters and the colour gradient)
            pushViewport(viewport(layout = grid.layout(2, 1, widths = c(1, 1), heights = c(0.1, 0.9))))
            pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
            pushViewport(viewport(layout = grid.layout(1, 1, widths = c(1, 1), heights = rep(1, 2))))
            grid.text(substitute(paste(k[max],"= NA, N = ", N, ", L = ", L, sep = ""), list(N=N, L=L)), x = 0.2, y = 0.65)
            grid.text(substitute(paste(delta, "= ", a, ", ",nu," = ", b, ", threshold = ", f, "%", sep = ""), list(a=d, b=v, f=threshold)), x = 0.2, y = 0.3)
            upViewport()
            pushViewport(viewport(layout = grid.layout(2, 2, widths = c(1, 1), heights = rep(1, 2))))

            pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
            pushViewport(plotViewport(margins = c(1.5, 1, 1, 1)))
            grid.points(seq(0, 1, length =  nrow(lists)), rep(1, nrow(lists)), pch = 15, gp = gpar(col = red.white.green))
            pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))	
            grid.text("close", x = 0, y = 0.5, gp = gpar(fontsize = 8, col = "black"))
            grid.text("RANK POSITION", x = 0.5, y = 0.5, gp = gpar(fontsize = 7, col = "black"))
            grid.text("distant", x = 1, y = 0.5, gp = gpar(fontsize = 8, col = "black"))
            upViewport()
            popViewport()
            popViewport()
            upViewport()
            upViewport()
                                        #set the layout 
            pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))	
            pushViewport(viewport(layout = grid.layout(1, 2, widths = c(1), heights = c(1,1))))
            grid.text(substitute(paste("On selected ", delta ,", no overlap found", sep="")))
            upViewport()
            upViewport()
            upViewport()  } #end for is.null 
}#end for aggmap


.isCalculatable <- function(lists, d, L, v) {
  res.j0.temp <- j0.multi(lists, d, v)
  res.temp <- as.matrix(res.j0.temp$L)
  
  temp2 = c()
                                        #calculate the reference lists
  for (i in c(1:L)) {
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
  
                                        #calculate the truncated lists for each block (for each reference list)
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
    
                                        #get the truncation points for the truncated lists in the current block
    for(currennt.truncated.list in ranked.list.names[2:length(ranked.list.names)]) {
      truncation.points <- c(truncation.points, as.numeric(res.temp[res.temp[,1] == first.list.name & res.temp[,2] == currennt.truncated.list,4]))
    }
    
    lists.to.remove <- c(lists.to.remove, first.list.name)
  }
  
                                        #check if a truncation point is missing - if true, then the delta-value cannot be applied on the lists (aggmap is not drawn)
  if(NA %in% truncation.points) {
    return(FALSE)
  }else {
    return(TRUE)
  }
}

