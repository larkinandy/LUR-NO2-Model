############ cross_validation.R ############
# Author: Andrew Larkin
# Developed for Perry Hystad, Oregon State University
# Date created: January 12, 2016
# This script performs a % leave out cross validation on the lassoo regression variable selection process for
# the global NO2 land use regression model for n number of times.  The user specific the number of bootstrap
# samples to run, and the remaining variables are hard coded into the script.  It is not recommended to 
# modify variables in the script as it is designed only to work with the specific NO2 dataset provided.

####### load required packages #########
library(glmnet) # lasso regression



######################## helper functions #####################


# for variables with protective effects, inverse the sign of the variable value
# to allow the lasso regression algorithm to force the direction of the coefficient
# to correspond to protective effects (reduction in NO2 concentrations)
# INPUTS:
#   inData (dataframe) - dataset containing all predictor variables
# OUPUTS:
#   inData (dataframe) - same dataset as supplied to the model, but with the 
#                        direction of the protective variables inversed
posCoeffMatrix <- function(inData) {
  
  tempNames <- names(inData)  # get the names of the variables
  endLength <- length(tempNames) # get the number of variables
  switchList <- c("tr","ND","wa","us","oe") # list which variables are protecctive
  
  # for each variable in the dataset
  for(i in 1:endLength) {
    predType <- substr(tempNames[i],1,2) # get the type of variable (major roads, NDVI, etc)
    if(predType %in% switchList) { # if the variable is in the list of protective variables, inverse the sign of the value
      inData[,i] <- inData[,i]*-1
    }                       
    
  }
  return(inData)
} # end of posCoeffMatrix




# calculate the root mean square error, mean abs error, r-square, adjusted r-square, bias, and abs bias for all regions
# INPUTS:
#   inResids (float vector) - residuals from the model prediction (residuals = monitor - prediction)
#   zoneVals (int vector) - the zone region corresponding to the residuals and air monitor measurements)
#   p (int) - number of variables in the lasso model
# OUTPUTS:
#   returnData (dataframe) - data frame containing the model evaluation summary statistics for each region
calcAllRMSE <- function(inResids,zoneVals,inMonitor,p) {
 
  uniqueZones <- unique(zoneVals)[order(unique(zoneVals))] # identify unique zones that need to be processed, and order
  mse <- ase <- rsq <- adjRsq <- percBias <- bias <- rep(0,length(uniqueZones)+1)
  
  # calculate summarry statistics for the entire dataset
  mse[7] <- sqrt(mean(inResids^2)) 
  ase[7] <- mean(abs(inResids))
  sumSqErr <- sum(inResids^2)
  sumTot <- sum((inMonitor - mean(inMonitor))^2)
  n <- length(inResids)
  rsq[7] <- 1 - (sumSqErr/sumTot)
  adjRsq[7] <- 1 - (((1-rsq[7])*(n-1))/(n-p-1))
  percBias[7] <- (100/length(inResids))*(sum(abs(inResids)/inMonitor))
  bias[7] <- (-100/length(inResids))*(sum(inResids/inMonitor))
  
  # for each region, calculate the summary statistics
  for(i in 1:length(uniqueZones)) {
    tempResid <- subset(inResids,zoneVals == uniqueZones[i])
    tempMonitor <- subset(inMonitor,zoneVals==uniqueZones[i])
    mse[i] <- sqrt(mean(tempResid^2))
    ase[i] <- mean(abs(tempResid))
    sumSqErr <- sum(tempResid^2)
    sumTot <- sum((tempMonitor - mean(tempMonitor))^2)
    n <- length(tempResid)
    rsq[i] <- 1 - (sumSqErr/sumTot)
    rsq[i] <- max(rsq[i],0)
    adjRsq[i] <- 1 - (((1-rsq[i])*(n-1))/(n-p-1))
    adjRsq[i] <- max(adjRsq[i],0)
    percBias[i] <- (100/length(tempResid))*(sum(abs(tempResid)/tempMonitor))
    bias[i] <- (-100/length(tempResid))*(sum(tempResid/tempMonitor))
  }
  
  # combine summary statistic vectors into a single dataframe to return
  returnData <- data.frame(mse,ase,rsq,adjRsq,bias,percBias)
  return(returnData)
} # end of calcAllRMSE


# create dummary intercept variables for each region
# INPUTS:
#    inData (dataframe) - contains the predictor variable dataset
#    zoneVals (int vector) - contains the zone that each row in the dataset corresponds to
# OUTPUTS:
#    inData (dataframe) - the input predictor variable dataset, with the region intercepts added to the end columns
#                         of the dataset
createRegionIntercepts <- function(inData,zoneVals) {

  inData$z1 <- (zoneVals == 1)*1
  inData$z3 <- (zoneVals == 3)*1
  inData$z4 <- (zoneVals == 4)*1
  inData$z6 <- (zoneVals == 6)*1
  inData$z7 <- (zoneVals == 7)*1
  inData$z9 <- (zoneVals == 9)*1
  return(inData)
} # end of createREgionalIntercepts




# using an input dataset, randomly partition a training and testing dataset for cross-validation
# INPUTS:
#    inData (dataframe) - input dataset with predictor variables and air monitor measurements
#    sampProp (float) - value ranging from 0 to 1, indicating the proportion of data that should be partitioned to the training dataset
#    zoneVals (int vector) - indicates which zone the corresponding row of the input dataset belongs to
# OutPUTS:
#    returnData (dataframe) - input dataset with an indicator variable of whether each row belongs to the train or test partition
createTrainAndTest <- function(inData,sampProp,zoneVals) {

  zoneVals1 <- subset(inData,zoneVals==1) # for the first zone, craete a train and test dataset
  smp_size <- floor(sampProp* nrow(zoneVals1)) # calculate the sample size for the training dataset
  train_ind <- sample(seq_len(nrow(zoneVals1)), size = smp_size) # randomly sample the entire dataset to make the training dataset
  train <- zoneVals1[train_ind, ]
  test <- zoneVals1[-train_ind, ]
  uniqueZoneVals <- c(3,4,6,7,9)
  
  # for each of the unique regions after region 1, repeat the random sampling process shown for zone 1 above
  for(i in 1:length(uniqueZoneVals)) {
    tempData <- subset(inData,zoneVals==uniqueZoneVals[i])
    smp_size <- floor(sampProp* nrow(tempData))
    train_ind <- sample(seq_len(nrow(tempData)), size = smp_size)
    train <- rbind(train,tempData[train_ind, ])
    test <- rbind(test,tempData[-train_ind, ])
  }
  
  # create an indicator variable for whether a given sample (row in the dataset is a train or test point)
  train$ind <- rep(0,nrow(train))
  test$ind <- rep(1,nrow(test))
  
  returnData <- rbind(train,test) # combine the train and test dataset and return the result
  return(returnData)
} # end of createTrainAndTest



# perform leave 10% out cross-validation numRep number of times.  Return the root mean square, mean abs square, r-square,
# adjusted r-square, bias, and abs bias
# INPUTS:
#    inData (dataframe) - input data frame containing both the predictor and air monitor variables
#    numReps (int) - number of cross-validation repititions to perform
# OUTPUTS:
#    returnData (dataframe) - summary statistics of the cross-validation, for each region
crossValidation <- function(inData,numReps) {
  
  rmse <- ase <- rsq <- adjRsq <- bias <- absBias <- rep(0,7)
  
  
  ###### exactMonitors ############
  exactMonitors <- subset(inData,exact==1)
  exactMonitors <- addZoneField(exactMonitors)
  exactMonitors <- subset(exactMonitors,numMeas >1 )
  exactMonitors <- subset(exactMonitors,minYr > 2003 )
  exactLat <- exactMonitors$latitude 
  exactLong <- exactMonitors$longitude
  exactMonitors$monitor <- exactMonitors$meanNO2
  
  # reduce variables in the dataset to only predictor variables
  drops <- c("latitude","longitude","minYr", "maxYr", "medYr", "minNO2", "meanNO2", "maxNO2", "stdDevNO2", "numMeas","dist2Coast",
             "exact","X","dist2Coast","decLat","decLong","FID","NAME","CONTINENT")
  exactMatrix<- exactMonitors[ , !(names(exactMonitors) %in% drops)]
  exactMatrix2 <- createRegionIntercepts(exactMatrix,exactMonitors$zone) # create regional dummary variables

  # for each cross-validation repitition
  for(i in 1:numReps) {
    cat(i) # print the repitition number to the screen
    cat("\n")
    trainInd <- createTrainAndTest(exactMatrix2,0.8,exactMonitors$zone) # create training and testing datasets
    
    # partition trainInd into training and test datasets based on the indicator variable ind ( 0 = 1, 1 = test)
    monitor <- trainInd$monitor 
    zone <- trainInd$zone
    zone <- subset(zone,trainInd$ind==1)
    trainInd <- trainInd[ , !(names(trainInd) %in% c("monitor","zone"))]
    trainSet <- subset(trainInd,ind == 0)
    trainMonitor <- subset(monitor,trainInd$ind == 0)
    testSet <- subset(trainInd,trainInd$ind == 1)
    testMonitor <- subset(monitor,trainInd$ind == 1)
    drops <- c("ind")
    trainSet <- trainSet[ , !(names(testSet) %in% c("ind"))]
    testSet <- testSet[ , !(names(testSet) %in% c("ind"))]
    
    # inverse the direction of the variable value for protective variables
    trainMat <- as.matrix(posCoeffMatrix(trainSet))
    testMat <- as.matrix(posCoeffMatrix(testSet))
    
    # perform lasso variable selection
    cvfit <-  glmnet::cv.glmnet(trainMat,trainMonitor,type.measure = "mse",standardize=TRUE,alpha = 1,lower.limit=0)
    coefRaw <- coef(cvfit, s = "lambda.1se")
    
    keeps <- reduceLassoModel(coefRaw,trainSet,3)
    
    trainSet2 <- trainSet[, (names(trainSet) %in% keeps)]
    exactMat2 <- as.matrix(posCoeffMatrix(trainSet2))
    
    
    cvfit <- glmnet::cv.glmnet(exactMat2,trainMonitor,type.measure = "mse",standardize=TRUE,alpha = 1, lower.limit=0)
    coefRaw<- coef(cvfit, s = 0)
    
    testSet2 <- testSet[,(names(testSet) %in% keeps)]
    testMat <- as.matrix(posCoeffMatrix(testSet2))
    testMat2 <- cbind(rep(1,nrow(testMat)),testMat)
    
    
    # create predictions for the test dataset based on the variables selected by the training dataset
    pred <- as.vector(coefRaw[1:length(coefRaw)]%*%t(testMat2))
    residuals <- testMonitor-pred
    
    # evaluate the results of the test dataset cross-validation
    evalData <- calcAllRMSE(residuals,zone,testMonitor,ncol(testMat2)-6)
    rmse <- rmse + evalData[,1]
    ase <- ase + evalData[,2]
    rsq <- rsq + evalData[,3]
    adjRsq <- adjRsq + evalData[,4]
    bias <- bias + evalData[,5]
    absBias <- absBias + evalData[,6]
    
  }
  
  returnData <- data.frame(rmse,ase,rsq,adjRsq,bias,absBias) # combine evaluation statistics into a dataframe to return as output
  return(returnData)
  
} # end of crossValidation


############ main function ############


setwd("Insert Working Directory path here")
rawData <- read.csv("insert csv filename within working directory here")

numReps <- 1000
sumRMSE <- crossValidation(rawData,numReps)
avgRMSE <- sumRMSE/numReps







##### end of cross_validation.R ######
