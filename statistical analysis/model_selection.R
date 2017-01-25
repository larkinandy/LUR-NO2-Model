############ model_selection.R ############
# Author: Andrew Larkin
# Developed for Perry Hystad, Oregon State University
# Date created: January 23, 2016
# This script performs lasso varaible selection and incremental varaible buffer reduction for the NO2 LUR model.
# RMSE, AME, R2, Adj. R2, MB, and MAB are calculated for the final model.  The model coefficients, and minimum 
# p-value and percent variance explained for each variable and each region are also calculated.

####### load required packages #########
library(glmnet) # lasso regression

######################## helper functions #####################



# create a matrix in which the sign of protective variables "tr, ND, wa, us, and oe" are flipped.  
# flipping the sign of protective variables allows the lasso regression to restrict coefficients to 
# positive coefficients only: that is, a positive coefficient of an inverted protective value is 
# equivalent to a negative value before the sign of the variable was flipped
# INPUTS:
#    inData (dataframe) - matrix containing dataset 
# OuTPUTS:
#    inData (dataframe) - same matrix as input data, but with the signs of the protective
#                         variables flipped
posCoeffMatrix <- function(inData) {
  tempNames <- names(inData)
  endLength <- length(tempNames)
  switchList <- c("tr","ND","wa","us","oe") # list of two characters that indicate protective variables
  # for each variable in the dataset, check if the variable is in the list of protective variables.
  # if the variable is protective, multiply the value by 1
  for(i in 1:endLength) {
    predType <- substr(tempNames[i],1,2)
    
    if(predType %in% switchList) {
      inData[,i] <- inData[,i]*-1
    }                     
  }
  return(inData)
  
} # end of posCoefMatrix



# create a multipanel image in ggplot2.  thanks to winston@stdout.org for providing this function.
# function was downloaded from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL, save) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if(save!=FALSE) {
    ppi <- 300
    png(save, width=10*ppi, height=12*ppi,res=ppi)
  }
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
        print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                        layout.pos.col = matchidx$col))
    }
  }
  
  if(save!=FALSE) {
    dev.off()
  }
} # end of multiplot







# add a categorical variable indicating which continent that monitor data is from
# INPUTS:
#    inputData (dataframe) - matrix containing dataset 
# OuTPUTS:
#    outputList (dataframe) - same matrix as input data, but with an added region
#                             variable
addZoneField <- function(InputData) {
  
  
  # create a list of all continents, along with a vector indicating which number
  # corresponds to which continent
  continentList <- c("North America","South America","Europe", "Africa","Asia")
  zoneList <- c(1,3,4,6,7)

  # combine Oceania and Australia together to create a single oceania continental region
  # use OCeania as the seed for the output dataset
  outputList <- subset(InputData, CONTINENT == "Oceania" | CONTINENT == "Australia")
  tempZone <- rep(9,length(outputList[,1]))
  outputList$zone <- tempZone
  
  # for each continent except oceania, add the subset to the output dataset along with 
  # the corresponding continental region value
  for(index in 1:length(continentList)) {
    tempData <- subset(InputData,CONTINENT == continentList[index])
    tempData$zone <- zoneList[index]
    outputList <- rbind(outputList,tempData)
  }
  
  return(outputList)
  
} # end of addZoneField




# identify the buffer sizes included in the variales in the input dataset.  This is done
# by removing the first two characters from each variable and converting the remaining 
# characters from characters to an integer value
# INPUTS:
#    inputData (dataframe) - matrix containing variables with buffer distances of interest 
# OuTPUTS:
#    buffDist (integer array) - array containing buffer distances, in ascending order
getBuffDistVec <- function(inputData) {
  buffDist <- rep(100,length(inputData)) # array that will contain output data
  # for each variable, extract the buffer distance from the variable name and conver to an
  # integer
  for(j in 1:length(inputData)) { 
    endP <- nchar(inputData[j]) -1 
    buffDist[j] <- as.numeric(substr(inputData[j],3,endP))
  } 
  buffDist <- buffDist[order(buffDist)] #order buffer distances in ascending order
  return(buffDist)
} # end of getBuffDistVec



# reduce buffers to only those that are more than x fold distance apart from one another
# INPUTS:
#    inputData (dataframe) - matrix containing variables with buffer distances of interest 
# OuTPUTS:
#    buffDist (integer array) - array containing buffer distances, in ascending order
reduceBuffList <- function(inputData,fold=3) {
  index = 1
  while(index<length(inputData)) {
    finishedVarCompare = FALSE
    while(finishedVarCompare == FALSE & index <length(inputData)) {
      if(inputData[index+1]/inputData[index] <= fold) {
        inputData <- inputData[-(index+1)]
      }
      else {
        finishedVarCompare=TRUE
      }
    }
    index = index +1
  }
  return(inputData)
} 



# reduce incremental variables within x fold values of the smallest vriable size
# INPUTS:
#    inCoef (float array) 
# OuTPUTS:
#    buffDist (integer array) - array containing buffer distances, in ascending order
reduceLassoModel <- function(inCoef,inPred,fold=5) {

  # create a vector of the first two characters for all variables
  bufferTypes <- c("fi","im","mj","mi","po","ND","tr","pl","uj","ui") 

  a <- which(inCoef > 0) # identify which variables were selected by lasso regression
  b <- a[2:length(a)]-1 # remove the intercept from the model 
  subNames <- names(inPred)[b] # get the names of the variables selected by lasso regression
  finalList <- c()
  index <- 0
  
  # for each type of varaible, remove variables that are within 3 fold of a smaller variable
  for(index in 1:length(bufferTypes)) {
    tempData <- subNames[substr(subNames,1,2) %in% bufferTypes[index]] 
    
    # get the the distances for all buffers of the selected variable type
    if(length(tempData)>0) {
      buffList <- getBuffDistVec(tempData)
      reduced <- reduceBuffList(buffList,fold)
      m <- "m"
      reduced <- paste(bufferTypes[index],reduced,m,sep="")
      finalList <- c(finalList,reduced)
    }
  }
  
  otherVars <- subNames[substr(subNames,1,2) %in% bufferTypes == FALSE]
  finalList <- c(finalList,otherVars)
  return(finalList)
} # end of reduceLassoModel



# get the summary statistics for each region.  Summary statistics include RMSE, AME, R2, Adj R2, MB, and MAB
# INPUTS:
#    inData (float array) - array of residuals
#    inZone (int array) - array of zone labels for other input variables
#    inMonitor (float array) - monitor measurements
#    p (integer) - number of variables in the model that corresponds to the residuals
# OuTPUTS:
#    returnData (dataframe) - dataframe containing summary statistics.  Each row corresponds to a 
#                             unique continental region
residualsByZone <- function(inData,inZone,inMonitor,p) {
  
  uniqueZones <- unique(inZone)[order(unique(inZone))] # get the unique zones in the dataset
  RMSE <- rsq <- adjRsq <- ase <- percBias <- bias <- rep(0,length(uniqueZones)+1)
  
  # for each unique zone, subset the residuals and monitors that correspond to the zone.  Then
  # calculate the summary statistics for that zone
  for(i in 1:length(uniqueZones)) {
    tempData <- subset(inData,inZone==uniqueZones[i])
    tempMonitor <- subset(inMonitor,inZone == uniqueZones[i])
    RMSE[i] <- sqrt(mean(tempData^2))
    sumSqErr <- sum(tempData^2)
    sumTot <- sum((tempMonitor - mean(tempMonitor))^2)
    n <- length(tempData)
    rsq[i] <- 1 - (sumSqErr/sumTot)
    adjRsq[i] <- 1 - (((1-rsq[i])*(n-1))/(n-p-1))
    ase[i] <- mean(abs(tempData))
    percBias[i] <- (100/length(tempData))*(sum(abs(tempData)/tempMonitor))
    bias[i] <- (-100/length(tempData))*(sum(tempData/tempMonitor))
  }
  
  # calculate summary statistics for the entire global dataset
  RMSE[length(RMSE)] <- sqrt(mean(inData^2))
  meanMonitor <- mean(inMonitor)
  sumSqErr <- sum(inData^2)
  sumTot <- sum((inMonitor - meanMonitor)^2)
  rsq[length(rsq)] <- 1 - (sumSqErr/sumTot)
  n <- length(inMonitor)
  adjRsq[length(rsq)] <- 1 - (((1-rsq[length(rsq)])*(n-1))/(n-p-1))
  ase[i+1] <- mean(abs(inData))
  percBias[i+1] <- (100/length(inData))*(sum(abs(inData)/inMonitor))
  bias[i+1] <- (-100/length(inData))*(sum(inData/inMonitor))
  returnData <- data.frame(RMSE,ase,rsq,adjRsq,bias,percBias,c(uniqueZones,0))
  return(returnData)
} # end of residualsByZone
    





# create a multipanel plot of the residuals for each continental region.  
# INPUTS:
#    predictions (float array) - array containing the model predictions
#    zoneVals (int array) - array containing zone values for predictions and monitors
#    monitor (float array) - array containing the monitor values
#    outputFilename (string) - output file name and filepath
#    plotType (string) - type of plot to make
makeResidPlots <- function(predictions,zoneVals,monitor,outputFileName) {
  uniqueZones <- unique(zoneVals)[order(unique(zoneVals))] # get unique zone values
  
  # if the only unique zone values is 9999, create a single global images
  if(uniqueZones==-9999) { 
    
    # define graph variables
    zoneNames <- c("Global")
    yLab = "monitor aveage (ppb)"
    xLab <- "predicted (ppb)"
    maxVal <- max(max(predictions),max(monitor))
    minVal <- min(min(predictions),min(monitor))
    predictions2 <- data.frame(predictions,monitor)
    
    # create the graph
    ggplot(data = predictions2, aes(x=predictions, y=monitor)) + geom_point() + 
    ggtitle(zoneNames[1]) + labs(x=xLab,y=yLab) + ylim(minVal,maxVal) + xlim(0,maxVal) + 
    geom_smooth(method = "lm", se = FALSE) + geom_abline(intercept=0,slope=1)
    #dev.off()
  }
  
  # create subset plots for each continental region
  else {
  zoneNames <- c("North America","South America","Europe","Africa","Asia","Oceania")
  par(mfrow = c(3, 2))
  yLab = "monitor aveage (ppb)"
  
  # for each continental region, create sub plots
  for(i in 1:length(uniqueZones)) {
    tempResid <- subset(predictions,zoneVals==uniqueZones[i])
    tempMonitor <- subset(monitor,zoneVals == uniqueZones[i])
    xLab <- "predicted (ppb)"
    maxVal <- max(max(tempResid),max(tempMonitor))
    minVal <- min(min(tempResid),min(tempMonitor))
    tempData <- data.frame(tempResid,tempMonitor)
    tempPlot <- ggplot(data = tempData, aes(x=tempResid, y=tempMonitor)) + geom_point() + 
      ggtitle(zoneNames[i]) + labs(x=xLab,y=yLab) + ylim(minVal,maxVal) + xlim(0,maxVal) + 
      geom_smooth(method = "lm", se = FALSE) + geom_abline(intercept=0,slope=1)
    assign(paste("p",as.character(i),sep=""), tempPlot)
  }
    #plot(tempMonitor,tempResid,xlim=c(minVal,maxVal),ylim=c(minVal,maxVal),title=)
  multiplot(p1,p2,p3,p4,p5,p6, cols = 2,save = outputFileName)
  }
} # end of makeResidPlots



# create country dummy variables
# INPUTS:
#    inData (dataframe) - matrix containing input predictor matrix
#    zoneVals (int array) - array indicating continental region for each row
# OUTPUTS:
#    inData (dataframe) - same as the inData but with country intercept variables 
#                         added
createCountryIntercepts <- function(inData,zoneVals) {
  inData$z1 <- (zoneVals == 1)*1
  inData$z3 <- (zoneVals == 3)*1
  inData$z4 <- (zoneVals == 4)*1
  inData$z6 <- (zoneVals == 6)*1
  inData$z7 <- (zoneVals == 7)*1
  inData$z9 <- (zoneVals == 9)*1
  return(inData)
} # end of createCountryIntercepts




# calculate and graph the partial R2 values for all variables and continental regions.
# INPUTS:
#    inData (dataframe) - data matrix containing predictor variables
#    inMonitor (float array) - array containing monitor measurements
#    outFile (string) - filepath and name of the output .eps file
graphPartialR2 <- function(inData, inMonitor, outFile) {
  uniqueZoneVals <- c(0,1,3,4,6,7,9) # list of all of the unique continental regions
  
  # for each continental region, subest the predictor matrix and monitors to the region of interest
  for(zoneVal in 1:length(uniqueZoneVals)) {
    
    if(uniqueZoneVals[zoneVal] == 0) {
      tempMat <- as.matrix(inData)
      tempMonitor <- inMonitor
    }
    
    else {
      tempMat <- subset(as.matrix(inData),exactMonitors$zone == uniqueZoneVals[zoneVal])
      tempMonitor <- subset(inMonitor,exactMonitors$zone == uniqueZoneVals[zoneVal])
      cat(zoneVal)
    }
    
    partialR2 <- rep(0,length(inData))
    lmTotal <- lm(tempMonitor~tempMat)
    
    # claculate partial R2 for all variables in the dataset
    ssrTot <- sum(anova(lmTotal)$"Sum Sq"[1:2])
    sseTot <- anova(lmTotal)$"Sum Sq"[2]
    for(i in 1:length(inData)) {
      tempRemove <- names(inData)[i]
      tempData <- tempMat[ , !names(inData) %in% tempRemove]
      tempLm <- lm(tempMonitor~tempData)
      tempSSR <- anova(tempLm)$"Sum Sq"[1]
      tempSSE <- anova(tempLm)$"Sum Sq"[2]
      partialR2[i] <- (tempSSE - sseTot)/tempSSE
      
      if(zoneVal == 1) {
        cat(partialR2[i])
        cat("\n")
      }
    }
    
    if(uniqueZoneVals[zoneVal]==0) {
      valMat <- round(partialR2*100,2)
      
    }
    else {
      valMat <- cbind(valMat,round(partialR2*100,2))
      cat(length(valMat[,1]))
      cat("\n")
    }
  }
  
  library(gplots)
  library(RColorBrewer)
  
  valMat <- data.frame(valMat[1:15,])
  names(valMat) <-c("Global","N America","S America","Europe","Africa","Asia","Oceania")
  rownames(valMat) <- names(inData)[1:15]
  
  
  # define parameters for creating heatmap
  
  "#FF2400"
  my_palette <- colorRampPalette(c("#E62020", "yellow", "green"))(n = 299)
  
  
  # (optional) defines the color breaks manually for a "skewed" color transition
  col_breaks = c(seq(0,0.49,length=100),   # for red
                 seq(0.5,1.49,length=100),            # for yellow
                 seq(1.5,2.5,length=100))              # for green
  
  row_distance = dist(valMat, method = "manhattan")
  row_cluster = hclust(row_distance, method = "ward.D")
  col_distance = dist(t(valMat), method = "manhattan")
  col_cluster = hclust(col_distance, method = "ward.D")
  
  setEPS()
  postscript(outFile)
  
  heatmap.2(as.matrix(valMat),
            cellnote = valMat,  # same data set for cell labels
            main = "Partial R2 2 Fold", # heat map title
            notecol = "black",      # change font color of cell labels to black#
            density.info = "none",  # turns off density plot inside color legend
            trace = "none",         # turns off trace lines inside the heat map
            margins = c(12,9),     # widens margins around plot
            col = my_palette,   
            breaks=col_breaks,# use on color palette defined earlier
            Rowv = as.dendrogram(row_cluster), # apply default clustering method
            Colv = as.dendrogram(col_cluster), # apply default clustering method
            dendrogram = "none")
  
  dev.off()
} # end of graphPartialR2






# calculate the minimum p-value of each variable for all continental and global regions
# INPUTS:
#    inMatrix (dataframe) - data matrix containing predictor variables
#    inMonitor (float array) - array containing monitor measurements
#    inZones (int array) - array indicating which zones each row of data belongs to
# OUTPUTS:
#    minPVals (int array) - array containing minimum p-values for each variable, in the same
#                           order as the predictor variables in the inMatrix
calcMinPValue <- function(inMatrix,inMonitor,inZones) {
  uniqueZones <- c(1,3,4,6,7,9) # identify unique continental regions
  minPVals <- rep(100,14) # create vector for p vals
  
  # for each continental region, subset the predictor matrix to the continental region, and 
  # calculate the pvals.  If the p vals are less than current min pvals, then update the pvals
  for(i in 1:length(uniqueZones)) {
    tempData <- subset(inMatrix,inZones==uniqueZones[i])
    tempMonitors <- subset(inMonitor,inZones==uniqueZones[i])
    tempModel <- summary(lm(tempMonitors ~ tempData))
    pVals <- tempModel$coefficients[,4]
    for(j in 1:length(pVals)) {
      if(pVals[j]<minPVals[j]) {
        minPVals[j] <- pVals[j]
      }
    }
  }
  return(minPVals)
} # end of calcMinPValue



# calculate the IQR for all variables in a matrix
# INPUTS:
#   inMatrix (dataframe) - variables for which IQR should be calculated
# OUTPUTS:
#   IQRVals (float array) - array of IQR values, in the same order as the 
#                           variables in inMatrix
calcIQR <- function(inMatrix) {
  IQRVals <- rep(0,length(inMatrix[1,]))
  for(i in 1:length(inMatrix[1,])) {
    IQRVals[i] <- IQR(inMatrix[,i])
  }
  return(IQRVals)
} # end of calcIQR


################# main script #################

library(ggplot2)
library(glmnet)

setwd("insert working directory here")



rawData <- read.csv("insert csv file with data here")

# setup data for processing
exactMonitors <- addZoneField(rawData)
exactMonitors <- subset(exactMonitors,numMeas >1 )
exactMonitors <- subset(exactMonitors,minYr > 2003 )
exactLat <- exactMonitors$latitude
exactLong <- exactMonitors$longitude
monitor <-  exactMonitors$meanNO2
drops <- c("zone","X","FID","decLat","decLong","exact","latitude","longitude","minYr", "maxYr", "medYr", "minNO2", "meanNO2", "maxNO2", "stdDevNO2", "numMeas","FID","NAME","CONTINENT")
exactMatrix <- exactMonitors[ , !(names(exactMonitors) %in% drops)]



######## run analysis on final model #########


keeps <- c("ND200m","ND1200m","im7000m","im1500m","wa50000m","mj100m","mj2500m","po3500m","tr1500m","satNO2","z1","z3","z4","z6","z7","z9")


exactMatrix2 <- createCountryIntercepts(exactMatrix,exactMonitors$zone) # create model with intercepts
exactMatrix2 <- exactMatrix2[, (names(exactMatrix2) %in% keeps)] # reduce variables to final model structure
tempMat <- as.matrix(posCoeffMatrix(exactMatrix2)) # reverse direction of protective variabless

cvfit <- glmnet::cv.glmnet(tempMat,monitor,type.measure = "mse",standardize=TRUE,alpha = 1,lower.limit=0) # perform lasso regression
coefRaw <- coef(cvfit,0)
exactMat3 <- cbind(rep(1,length(tempMat[,1])),tempMat)

# create predictions baesd on the regression
pred <- as.vector(coefRaw[1:length(coefRaw)]%*%t(exactMat3))

# calculate residuals and export to csv
residuals <- monitor-pred
residualsByZone(pred,exactMonitors$zone,monitor,13)
residualData <- cbind(exactLat,exactLong,residuals)
write.csv(residualData,"C:/users/user/desktop/residuals.csv")


# graph predicted variables vs. monitor observations
a<- data.frame(pred,monitor)
ggplot(data= a,aes(x=pred,y=monitor)) + geom_point() + ylim(0,60) + xlim(0,60) + 
  geom_abline(intercept=0,slope=1)+
ggtitle("a") + labs(x="Predicted (ppb)",y="Monitor average (ppb)") + ylim(0,60) + xlim(0,60) + 
  geom_smooth(method = "lm", se = FALSE) + geom_abline(intercept=0,slope=1)

# make residual plots for each continental region
makeResidPlots(pred,exactMonitors$zone,monitor,"D:/residPlotsCountryIntercept2.png")

# calculate the minimum p value 
calcMinPValue(tempMat,monitor,exactMonitors$zone)

######## create intercept models ############

exactMatrix2 <- createCountryIntercepts(exactMatrix,exactMonitors$zone) # add intercept variables to model structure
tempMat <- as.matrix(posCoeffMatrix(exactMatrix2)) # inverse the direction of protective variables

# create a lasso regression model
cvfit <- glmnet::cv.glmnet(tempMat,monitor,type.measure = "mse",standardize=TRUE,alpha = 1,lower.limit=0)

# identify varaibles that should be kept in the second step of the model building
coefRaw <- coef(cvfit, s = "lambda.1se")
keeps <- reduceLassoModel(coefRaw,exactMatrix,3)

# create model predictions and calculate residuals
exactMat3 <- cbind(rep(1,length(tempMat[,1])),tempMat)
pred <- as.vector(coefRaw[1:length(coefRaw)]%*%t(exactMat3))
residuals <- monitor-pred

# reduce model structure based variables that were identified to keep
tempData2 <- exactMatrix[, (names(exactMatrix) %in% keeps)]
tempData2 <- createCountryIntercepts(tempData2,exactMonitors$zone)
exactMat2 <- as.matrix(posCoeffMatrix(tempData2))

# refit the lasso regression model based on selected variables
cvfit <- glmnet::cv.glmnet(exactMat2,monitor,type.measure = "mse",standardize=TRUE,alpha = 1, lower.limit=0)
coefRaw<- coef(cvfit, s = 0)
exactMat3 <- cbind(rep(1,length(exactMat2[,1])),exactMat2)

# create predictions based on the updated model structure and calculate residuals
pred <- as.vector(coefRaw[1:length(coefRaw)]%*%t(exactMat3))
residuals <- monitor-pred

# calculate summary statistics for the global dataset and each continental region
residualsByZone(residuals,exactMonitors$zone,monitor,ncol(exactMat2))



########### end of ModelSelection.R ##########