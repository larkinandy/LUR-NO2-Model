############ sensitivity_analysis.R ############
# Author: Andrew Larkin
# Developed for Perry Hystad, Oregon State University
# Date created: January 24, 2016
# This script performs regional and residual model sensitivity analyses for the global NO2 LUR model.  
# Regional models are created based on continental subsets of the global dataset.  Residual models
# are based on creating LUR model to fit residuals for each continental region.    


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

    




calcSummaryStats <- function(inData,inZone,inMonitor,statType='mse') {
  zoneLabels <- c("North America","South America", "Europe", "Africa", "Asia", "Oceania","Global")
  if(statType == "mse") {
    allSigResid <- calcMSE(inData,inZone,inMonitor,lowerLim=-99999,reduceType=FALSE)
    valSelectResid <- calcMSE(inData,inZone,inMonitor,lowerLim=-99999,reduceType = TRUE)
    forcedResid <- calcMSE(inData,inZone,inMonitor,lowerLim=0,reduceType=TRUE)
    returnData <- data.frame(allSigResid,valSelectResid,forcedResid)
  }
  else if(statType == "r2") {
    allSigR2 <- calcR2(inData,inZone,inMonitor,lowerLim=-99999,reduceType=FALSE)
    valSelectR2 <- calcR2(inData,inZone,inMonitor,lowerLim=-99999,reduceType = TRUE)
    forcedR2 <- calcR2(inData,inZone,inMonitor,lowerLim=0,reduceType=TRUE)
    returnData <- data.frame(allSigR2,valSelectR2,forcedR2)
  }
    return(returnData)
}

calcRegional <- function(inData, zoneVals, monitor,zoneNum=0,lowerLim=0,reduce=FALSE,returnType="cvfit") {
  if(zoneNum>0) {
  tempData <- subset(inData, zoneVals ==zoneNum)
  tempMonitors <- subset(monitor,zoneVals==zoneNum)
  }
  else {
    tempData <- inData
    tempMonitors <- monitor
  }
  tempMat <- as.matrix(posCoeffMatrix(tempData))
  cvfit <- glmnet::cv.glmnet(tempMat,tempMonitors,type.measure = "mse",standardize=TRUE,alpha = 1,lower.limit=lowerLim)
  if(reduce==TRUE) {
    coefRaw <- coef(cvfit, s = "lambda.1se")
    cat(zoneNum)
    keeps <- reduceLassoModel(coefRaw,tempData,5)
    cat(keeps)
    tempData2 <- tempData[, (names(tempData) %in% keeps)]
    
    exactMat2 <- as.matrix(posCoeffMatrix(tempData2))
    cvfit <- glmnet::cv.glmnet(exactMat2,tempMonitors,type.measure = "mse",standardize=TRUE,alpha = 1, lower.limit=lowerLim)
  }
  else {
    tempData2 <- tempData
    exactMat2 <- as.matrix(posCoeffMatrix(tempData2))
  }
  if(returnType=="cvfit") {
    return(cvfit)
  }
  else if(returnType == "residual") {
    coefRaw<- coef(cvfit, s = "lambda.1se")
    exactMat3 <- cbind(rep(1,length(exactMat2[,1])),exactMat2)
    pred <- as.vector(coefRaw[1:length(coefRaw)]%*%t(exactMat3))
    residuals <- tempMonitors-pred
    return(residuals)
  }
  else if(returnType == "mse") {
    coefRaw<- coef(cvfit, s = "lambda.1se")
    exactMat3 <- cbind(rep(1,length(exactMat2[,1])),exactMat2)
    pred <- as.vector(coefRaw[1:length(coefRaw)]%*%t(exactMat3))
    residuals <- tempMonitors-pred
    mse <- sqrt(mean(residuals^2))
    return(mse)
  }
  else {
    return(cvfit$glmnet.fit$dev.ratio[length(cvfit$glmnet.fit$dev.ratio)])
  }
    
}

createCountryIntercepts <- function(inData,zoneVals) {
  inData$z1 <- (zoneVals == 1)*1
  inData$z3 <- (zoneVals == 3)*1
  inData$z4 <- (zoneVals == 4)*1
  inData$z6 <- (zoneVals == 6)*1
  inData$z7 <- (zoneVals == 7)*1
  inData$z9 <- (zoneVals == 9)*1
  return(inData)
}





# create a residual model for a given continental region, and return either summary statistics or the 
# model residuals
# INPUTS:
#    inZone (int) - integer indicating the continental region to create residual models for
#    returnType (string) - indicates which information about the residual model to return.  Default is mse, which 
#                          returns summary statistics.  
calcResidualModel <- function(inZone,returnType = "mse") {
  
  keeps <- c("ND100m","ND1200m","im7000m","im1500m","wa50000m","mj100m","mj2500m",
             "uj2500m","po3500m","tr1500m","satNO2")
  
  # reduce dataset to final model variables
  exactMatrix2 <- exactMatrix[, (names(exactMatrix) %in% keeps)]
  exactMatrix2 <- createCountryIntercepts(exactMatrix2,exactMonitors$zone)
  
  # invert the direction of the protective variable signs
  exactMat2 <- as.matrix(posCoeffMatrix(exactMatrix2))
  # fit the model using lasso variable selection
  cvfit <- glmnet::cv.glmnet(exactMat2,monitor,type.measure = "mse",standardize=TRUE,alpha = 1, lower.limit=0)
  coefRaw<- coef(cvfit, s = "lambda.1se")
  
  # create model predictions
  exactMat3 <- cbind(rep(1,length(exactMat2[,1])),exactMat2)
  pred <- as.vector(coefRaw[1:length(coefRaw)]%*%t(exactMat3))
  
  # subset predictions to the region of interest and create residuals
  naPred <- subset(pred,exactMonitors$zone==inZone)
  residuals <- monitor-pred
  northAmResid <- subset(residuals,exactMonitors$zone==inZone)

  # create a new predictor matrix to fit the residual model
  tempData2 <- exactMatrix
  a <- posCoeffMatrix(tempData2)
  exactMat2 <- as.matrix(a)
  # subset predictor matrix to the regional area of interest and fit the residual model
  northAmModels <- subset(exactMat2,exactMonitors$zone == inZone)
  cvfit <- glmnet::cv.glmnet(northAmModels,northAmResid,type.measure = "mse",standardize=TRUE,alpha = 1, lower.limit=0)
  coefRaw <- coef(cvfit, s = "lambda.1se")
  
  # subset residual model to significant variables
  a <- which(coefRaw > 0)
  if(length(a)  <=1 ) {
    predSub <- rep(coefRaw[1],length(northAmResid))
  }
  else {
  keeps <- reduceLassoModel(coefRaw,tempData2,3) # resduce lasso model 
  tempData2 <- exactMatrix[, (names(exactMatrix) %in% keeps)]
  tempData2 <- subset(tempData2,exactMonitors$zone == inZone)
  
  # if the residual model only has one variable, use linear regression to determine the coefficient
  if(length(keeps)==1) {
    lmModel <- lm(northAmResid~ tempData2)
    predSub <- lmModel$fitted.values
    naMonitors <- subset(monitor,exactMonitors$zone==inZone)
    residuals <- naMonitors-predSub-naPred
    
    # calculate summary statistics 
    ase <- mean(abs(residuals))
    sumSqErr <- sum(residuals^2)
    sumTot <- sum((naMonitors - mean(naMonitors))^2)
    rsq <- 1 - (sumSqErr/sumTot)
    mse <- sqrt(mean(residuals^2))
    n <- length(residuals)
    p <- 1
    adjRsq <- 1 - (((1-rsq)*(n-1))/(n-p-1))
    absBias <- (100/length(residuals))*(sum(abs(residuals)/naMonitors))
    bias <- (-100/length(residuals))*(sum(residuals/naMonitors))
    if(returnType == "mse") {
      return(c(mse,ase,rsq,adjRsq,bias,absBias))
    }
    else if ( returnType == "residuals"){
      return(residuals)
    }
    else {
      VarTypes <- names(lmModel$coefficients)
      varVals <- lmModel$coefficients
      returnData <- data.frame(VarTypes,varVals)
      return(returnData)
    }
  }
  
  # refit the coefficients of the predictor matrix using the reduced variable dataset
  exactMat2 <- as.matrix(posCoeffMatrix(tempData2))
  cvfit <- glmnet::cv.glmnet(exactMat2,northAmResid,type.measure = "mse",standardize=TRUE,alpha = 1, lower.limit=0)
  coefRaw<- coef(cvfit, s = 0)
  
  # create predictions based on the new dataset
  exactMat2 <- cbind(rep(1,length(exactMat2[,1])),exactMat2)
  predSub <- as.vector(coefRaw[1:length(coefRaw)]%*%t(exactMat2))
  }
  
  # calculate residuals and model summary statistics
  naMonitors <- subset(monitor,exactMonitors$zone==inZone)
  residuals <- naMonitors-predSub-naPred
  sumSqErr <- sum(residuals^2)
  sumTot <- sum((naMonitors - mean(naMonitors))^2)
  rsq <- 1 - (sumSqErr/sumTot)
  mse <- sqrt(mean(residuals^2))
  n <- length(residuals)
  p <- max(1,length(which(coefRaw>0))-1)
  adjRsq <- 1 - (((1-rsq)*(n-1))/(n-p-1))
  ase <- mean(abs(residuals))
  absBias <- (100/length(residuals))*(sum(abs(residuals)/naMonitors))
  bias <- (-100/length(residuals))*(sum(residuals/naMonitors))
  if(returnType == "mse"){
  return(c(mse,ase,rsq,adjRsq,bias,absBias))
  }
  else if(returnType == "residuals") {
    return(residuals)
  }
  else {
    # return model structure
    a <- which(coefRaw > 0)
    if(length(a) > 1) {
    b <- a[2:length(a)] - 1
    varNames <- c("intercept",names(tempData2)[b])
    coeff <- coefRaw[which(coefRaw>0)]
    returnData <- data.frame(varNames,coeff)
    return(returnData)
    }
    else {
      return(coefRaw[1])
    }
  }
} # end of calcResidualModel




################# main script #################

library(ggplot2)
library(glmnet)

# set working directory
setwd("insert working directory here")

# read in raw data file
rawData <- read.csv("insert csv file with data here")



###### screen monitors using selecction criteria ############
exactMonitors <- subset(rawData,exact==1)
exactMonitors <- addZoneField(exactMonitors)
exactMonitors <- subset(exactMonitors,numMeas >1 )
exactMonitors <- subset(exactMonitors,minYr > 2003 )


exactLat <- exactMonitors$latitude
exactLong <- exactMonitors$longitude
monitor <- exactMonitors$meanNO2
drops <- c("exact","zone","decLat","decLong","latitude","longitude","minYr", "maxYr", "medYr", "minNO2", "meanNO2", "maxNO2", "stdDevNO2", "numMeas","FID","NAME","CONTINENT")
exactMatrix <- exactMonitors[ , !(names(exactMonitors) %in% drops)]



############# reidual models ###############
uniqueZoneVals <- c(1,3,4,6,7,9)
a <- calcResidualModel(1,"mse")
for(i in 2:length(uniqueZoneVals)) {
  a <- rbind(a,calcResidualModel(uniqueZoneVals[i],"mse"))
}


############## regional models ############

# set which zone is the value for
zoneVal <- 1

# subset the matrix and monitor to the zone of interest
subsetMatrix <- subset(exactMatrix,exactMonitors$zone==zoneVal)
subsetMonitors <- subset(monitor,exactMonitors$zone==zoneVal)

# invert coefficients of protective variables
tempMat <- as.matrix(posCoeffMatrix(subsetMatrix))
# fit the model using lasso regression
cvfit <- glmnet::cv.glmnet(tempMat,subsetMonitors,standardize=TRUE,alpha = 1,lower.limit=0)
coefRaw <- coef(cvfit, s = "lambda.1se")
# create model predictions and calculate residuals
exactMat2 <- cbind(rep(1,length(tempMat[,1])),tempMat)
pred <- as.vector(coefRaw[1:length(coefRaw)]%*%t(exactMat2))
residuals <- subsetMonitors-pred

keeps <- reduceLassoModel(coefRaw,subsetMatrix,3)

# reduce matrix to selected variables and invert protective variables
tempData2 <- subsetMatrix[, (names(subsetMatrix) %in% keeps)]
exactMat2 <- as.matrix(posCoeffMatrix(tempData2))

# repeat lasso regression with reduced matrix
cvfit <- glmnet::cv.glmnet(exactMat2,subsetMonitors,type.measure = "mse",standardize=TRUE,alpha = 1, lower.limit=0)

# create prediction matrix and calculate residuals
coefRaw<- coef(cvfit ,s=0)
exactMat3 <- cbind(rep(1,length(exactMat2[,1])),exactMat2)
pred <- as.vector(coefRaw[1:length(coefRaw)]%*%t(exactMat3))
residuals <- subsetMonitors-pred

# calculate summary statistics
residualsByZone(residuals,rep(zoneVal,length(residuals)),subsetMonitors,length(exactMat3[1,]))

