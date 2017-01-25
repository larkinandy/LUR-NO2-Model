######## get_wind_vectors_v2.py ##########
# download 6 hourly wind direction estimates for locations of interest from Google Earth Engine.  
# locations of interest are provided by latitude longitude coordinates in a .csv file.  Results are given
# as angular wind direction for all 6 hour measurements within the specific year in a .csv file
# Author: Andrew Larkin
# Developed for Perry Hystad, Oregon State University
# Date last modified: January 25, 2017

######### import necessary modules #########
import ee
import time
import datetime
import csv
import math
import os
import sys
import urllib2
ee.Initialize()


####### variables user might want to change ############
DisplayGUI = False
folder = os.path.dirname(sys.argv[0]) + "/"
# year of interest must be the last four characters in the name of the .csv file
defaultCSV = folder + "inputCSV/" + "NewDeliPoints2014.csv" 
outputFile = "outputfile"
bufferDistance = 50
BATCH_SIZE = 1000 # number of locations to run at one time.  Too many locations will cause GEE to timeout and thorw errors
surveyYear = int(defaultCSV[len(defaultCSV)-8:len(defaultCSV)-4])
startDate = str(surveyYear) +  '01-01'
endDate = str(surveyYear+1) + '01-01'
collectionName = 'NOAA/CFSV2/FOR6H' # unique identifier in GEE for the wind dataset of interest
tempCSVFolder = folder + "tempCSVs/"
selectVar = ["u-component_of_wind_height_above_ground","v-component_of_wind_height_above_ground"]
CSV_DICT = ["latitude","longitude"] # headers in input csv file that identify latitude, longitude variables
NEW_LINE = " \n "
defaultConstants = False
tempImage = ee.Image()



# read data from csv file into custom data structure
# INPUTS:
#   inputFile (string) - filepath to input .csv file containing coordinates
# OUTPUTS:
#   csvArray (float array) - custom data structure containing latitude,longitude coordinates
def readCSVFile(inputFile):
    csvArray = [[],[]]
    with open(inputFile) as f: # open inthe input file and move past the 6 lines of headers     
        reader = csv.DictReader(f, dialect=f) # identify the names of the columns (i.e. variables) in the csv file
        for row in reader: # read a row as {column1: value1, column2: value2,...}
            for (k,v) in row.items(): # go over each column name and value 
                if(k in CSV_DICT): # if the variable name is in the list of variables in linkWebDict (variables we want to keep)
                    csvArray[CSV_DICT.index(k)].append(float(v))
        del reader, row
        print("completed importing csv file")
        return(csvArray)     



# convert latitude longitude coordinates into point objects in Google Earth Engine
# INPUTS:
#   csvArray (float array) - custom python data structure containing lat/long coordinates
# OUTPUTS:
#   featureList (object array) - array of GEE point objects.  The GEE equivalent of a point shapefile
def makeGeomPoints(csvArray):
    i = 0
    featureList = []
    while i < len(csvArray[1]):
        featureList.append(ee.Feature(ee.Geometry.Point(csvArray[1][i],csvArray[0][i])))
        i+=1
    print ("completed creating geometric points")
    return(featureList)



# filter the wind database only to the time interval of interest
# INPUTS:
#   inStrat (string) - start date for the time interval of interest
#   inEnd (string) - end date for the time interval of interest
# OUTPUTS:
#   datedCollect (image collection) - reduced image collection 
def filterCatalogSet(inStart,inEnd):
    collection = ee.ImageCollection(collectionName)
    datedCollect = collection.filterDate(inStart,inEnd)
    return(datedCollect)



# extract values from a raster at points of interest
# INPUTS:
#   inputVals (point array) - array of GEE points
def bufferRunMap(inputVals):
    filteredCatalog = filterCatalogSet(startDate,endDate)
    catalogDates = filteredCatalog.toList(10000)
    timeStart = band.get('system:time_start') 
    a = band.reduceRegion(ee.Reducer.mean(),inputVals.geometry(),3)#,proj)
    inputVals2 = inputVals.set({str(timeStart.getInfo()):a})
    return(inputVals2)
  
  
  
# extract values from all wind rasters during the time interval of interest
# INPUTS:
#   inputPointDataset (array obj) - array of GEE points
#   inputRasterDataset (imageCollection) - raster catalog of wind rasters during the time interval of interest
# OUTPUTS:
#   currentDataset (array obj) array of GEE points, with wind direction appended to attribute table
#   timeStamps (float array) time values that correspond to the wind direction measurements appeneded to GEE point attribute table
def extractWindValues(inputPointDataset,inputRasterDataset,numRasters):
    print("extracting NDVI values from satellite imagery...")
    timeStamps = []
    print(numRasters.getInfo())
    for rasterNum in range(0,numRasters.getInfo()):
        if((rasterNum%2) ==0):
            time.sleep(0.25)
        tempImage = ee.Image(inputRasterDataset.get(rasterNum))
        if(rasterNum==0):
            currentDataset = inputPointDataset
        global band
        band = ee.Image(tempImage.select(selectVar))
        timeStart = band.get('system:time_start')
        timeStamps.append(str(timeStart.getInfo()))
        currentDataset = currentDataset.map(lambda f: bufferRunMap(f))
        del band    
        print("       completed " + str(rasterNum) + " out of " + str(numRasters.getInfo()) + " rasters")
    print("completed extracting NDVI values")
    return([currentDataset,timeStamps])



# import CSV file and transform into GEE point format
# INPUTS:
#   inputCSV (string) - filepath to the CSV file containing lat/long coords of interest
# OUTPUTS:
#   residPoints (obj array) - array of GEE point objects
#   csvArray (float obj) - custom python obj containing arrays of lat/long coords
def importCSVData(inputCSV):
    csvArray = readCSVFile(inputCSV)
    residPoints = makeGeomPoints(csvArray)    
    return([residPoints,csvArray])



# load and filter wind direction raster data for processing
# INPUTS:
#   inStartDate (string) - begining date of temporal coverage of interest
#   inEndDate (string) - end date of temporal coverage of interest
# OUTPUTS:
#   catalogDates (string array) - array of dates that correspond to dates of wind 
#                                 rasters in filtered catalog
#   numRasters (int) - number of rasters in the filtered catalog
def setupRasterData(inStartDate, inEndDate):
    filteredCatalog = filterCatalogSet(inStartDate, inEndDate)
    catalogDates = filteredCatalog.toList(3000)
    numRasters = ee.Number(filteredCatalog.size())    
    return([catalogDates,numRasters])



# extract data from GEE point objects and store extracted data in python lists
# INPUTS:
#   pointObject (obj array) - array of GEE point objects
#   timeStamps (float array) - array of time stamps, in epoch format
#   numPoints (int) - number of obervations in each point object to process
# OUPUTS:
#   pointStrings (array list) - list of arrays, in which each array contains extracted values for a single GEE point
def extractDataFromPoints(pointObject,timeStamps,numPoints):
    listVersion = pointObject.toList(100000)
    valList = []
    pointsDone = 0
    index = 0
    ee.data.setDeadline(0)
    params2={'fileFormat':'CSV','selectors':'#','filename':'testFile'}
    while(pointsDone < numPoints):
        pointStrings = []
        validDownload = False
        try:
            url = pointObject.getDownloadURL("csv",timeStamps)
            print(url)
            response = urllib2.urlopen(url)    
            cr = csv.reader(response)
            pointStrings = []
            i=0
            for row in cr:
                print(row)
                pointString = ""
                for col in row:
                    if(i==0):
                        pointString = pointString + col + ","
                    else:
                        commaIndex = col.find(',')
                        uVal = col[0:commaIndex]
                        uVal = str(uVal).replace("{u-component_of_wind_height_above_ground=","")
                        try:
                            uVal = float(uVal)
                        except:
                            uVal = float(-9999)
                        vVal = col[commaIndex+2:len(col)-1]
                        vVal = str(vVal).replace("v-component_of_wind_height_above_ground=","")
                        try:
                            vVal = float(vVal)
                        except:
                            vVal = float(-9999)
                        combVal = (450 - math.atan2(vVal, uVal)*57.29580) % 360
                        pointString = pointString + str(combVal) + ","
                pointString = pointString[0:len(pointString)-1] + NEW_LINE
                pointStrings.append(pointString)
                i+=1
            pointStrings = pointStrings[1:len(pointStrings)]
            pointsDone = i
        except Exception as e:
            print(str(e))            
    print("completed extracting values from geom points")
    return(pointStrings)
        
    
# create header for the output CSV file using timestamps
# INPUTS:
#   csvHeaderVals (int array) - array of timestamps in epoch time
# OUTPUTS:
#   timeString (string) - timestamps concatenated in a comma separated string
def getHeaderDates(csvHeaderVals):
    timeString = ""
    for epochTime in csvHeaderVals:    
        timeString = timeString + epochTime + ","
    timeString = timeString[0:len(timeString)-1]
    return(timeString)

# write CSV file containing data for a single month of observations
# INPUTS:
#   CSVHeader (string) - header for the csv file in string format
#   rowVals (string array) - array of strings, with each string corresponding to a single row in the CSV file
#   batchNumber (int) - number of the current batch 
#   folder (string) - filepath to the location where the CSV file will be written
def writeTempCSV(CSVHeader,rowVals,batchNumer,folder):
    i = 0
    tempFilename =folder + "tempCSVsBatch" + str(batchNumer) + ".csv"
    a = open(tempFilename, 'wb')
    a.write(CSVHeader)
    a.write(" \n ")
    for singleRow in rowVals:
        a.write(singleRow)
    a.close()
    
# monthly CSV files into annual file of wind direction
# INPUTS:
#   inputCSVFolder (string) - folder containing monthly CSV files
# OUTPUTS:
#   csvArray (matrix) - matrix containing latitude, longitude, wind direction and time stamps for an entire year
def mergeCSVFiles(inputCSVFolder):
    fileList = sorted(os.listdir(inputCSVFolder))
    
    tempDict = []
    csvArray = []
    for fileNum in range(0,len(os.listdir(inputCSVFolder))):
        fileName = 'tempCSVsBatch' + str(fileNum) + ".csv"
        print(fileName)
        with open(inputCSVFolder + fileName) as f: # open inthe input file and move past the 6 lines of headers     
            reader = csv.DictReader(f, dialect=f) # identify the names of the columns (i.e. variables) in the csv file
            i=0
            for row in reader: # read a row as {column1: value1, column2: value2,...}
                if(i==0):
                    for (k,v) in row.items(): # go over each column name and value 
                        tempDict.append(k.strip(' '))
                        tempArray = [k.strip(' ')]
                        csvArray.append(tempArray)
                for (k,v) in row.items():
                    csvArray[tempDict.index(k.strip(' '))].append(v)        
                i+=1
            del reader, row
        print("completed importing csv file")
    return(csvArray)     


# write data in custom array format to a csv file
# INPUTS:
#  inArray (float array) - float values to write to csv file
#  latLongArray (obj array) - data structure containing lat/long coords
#  inputFolder (string) - folderpath to write csv file to
#  batchNum (int) - number of current batch
def writeOutputArray(inArray, latLongArray,inputFolder,batchNum):
    resultsFilename =inputFolder + outputFile + str(batchNum) + ".csv"
    a = open(resultsFilename, 'wb')
    a.write("latitude,longitude,")
    inArray = sorted(inArray)
    for inCol in inArray:
        a.write(inCol[0])
        a.write(',')
    for i in range(1,len(inArray[0])-1):
        a.write(NEW_LINE)
        a.write(str(latLongArray[0][batchNum*BATCH_SIZE + i-1]))
        a.write(",")
        a.write(str(latLongArray[1][batchNum*BATCH_SIZE + i-1]))
        a.write(",")
        for j in range(0, len(inArray)):
            a.write(str(inArray[j][i]))
            if(j < len(inArray)):
                a.write(',')  
    a.close()    
    
    
# read batch CSV files and combine into an array object
# INPUTS:
#   inputFolder (string) - folder containing batch CSV files that need to be combined.
# OUTPUTS:
#   csvArray (string) - array data from individual CSV files combined into single data structure
def processMainCSVS(inputFolder):
    fileList = sorted(os.listdir(inputFolder))
    filepaths = []
    numFiles = 0
    for filename in fileList:
        if(filename[len(filename)-3:len(filename)] == "csv"):
            filepaths.append("outputfile" + str(numFiles) + ".csv")
            numFiles+=1
    
    tempDict = []
    csvArray = []
    j=0
    for fileName in filepaths:
        i=0
        if(j< BATCH_SIZE):
            print(fileName)
            with open(inputFolder + fileName) as f: # open inthe input file and move past the 6 lines of headers     
                reader = csv.DictReader(f, dialect=f) # identify the names of the columns (i.e. variables) in the csv file
                for row in reader: # read a row as {column1: value1, column2: value2,...}
                    if(i==0 and j==0):
                        for (k,v) in row.items(): # go over each column name and value 
                            tempDict.append(k.strip(' '))
                            tempArray = [k.strip(' ')]
                            csvArray.append(tempArray)
                    for (k,v) in row.items():
                        csvArray[tempDict.index(k.strip(' '))].append(v)        
                    i+=1
                del reader, row
            print(i)
            j+=1
            print("completed importing csv file")
    return(csvArray)     



# combine batch CSV files into final output csv
# INPUTS:
#   inArray (string array) - array containing data to write to CSV file
#   outputFilename (string) - filepath to write output file
def writeFinalOutput(inArray,outputFilename):
    resultsFilename = outputFilename
    a = open(resultsFilename, 'wb')
    for i in range(0,len(inArray[0])):
        for inCol in inArray:
            a.write(str(inCol[i]))
            a.write(',')        
        a.write(NEW_LINE)
    a.close()  



# main function
def main():
   
    csvFile = defaultCSV
    startDate = str(surveyYear) + "-01-01"
    endDate = str(surveyYear+1) + "-01-01"
    [entireCSVData,latLongArray] = importCSVData(csvFile)
    for batchNum in range(0,int(math.ceil(len(entireCSVData)/(BATCH_SIZE*1.0)))):
        batchFolder = tempCSVFolder + "b" + str(batchNum) + "/"
        if not os.path.exists(batchFolder):
            os.makedirs(batchFolder)
        for i in range(0,12):
            print(batchFolder + "tempCSVsBatch" + str(i) + ".csv")
            if(os.path.isfile(batchFolder + "tempCSVsBatch" + str(i) + ".csv")):
                print("batch " + str(i) + " has already been processed.  Moving on to next batch.")
            else:
                monthStart = i%12 + 1 
                monthEnd = (i+1)%12 + 1
                startDate = str(surveyYear) + "-" + str(monthStart) + "-" + "01"
                if(i==11):
                    endDate = str(surveyYear+1) + "-" + str(monthEnd) + "-" + "01"
                else:
                    endDate = str(surveyYear) + "-" + str(monthEnd) + "-" + "01"
                
                [rasterPartition,numRasters] = setupRasterData(startDate,endDate)
                batchSubsetCSV = entireCSVData[(batchNum*BATCH_SIZE):min((batchNum+1)*BATCH_SIZE,len(entireCSVData))]
                batchPoints = ee.FeatureCollection(batchSubsetCSV)
                [extractedPoints,timeStamps] = extractWindValues(batchPoints,rasterPartition,numRasters)
                pointStrings = extractDataFromPoints(extractedPoints,timeStamps,len(batchSubsetCSV))
                headerTimeString = getHeaderDates(timeStamps)
                writeTempCSV(headerTimeString,pointStrings,i,batchFolder)
        NDVIArray = mergeCSVFiles(batchFolder)
        writeOutputArray(NDVIArray,latLongArray,tempCSVFolder,batchNum)
    NDVIArray = processMainCSVS(tempCSVFolder)
    writeFinalOutput(NDVIArray, tempCSVFolder + "nycp42014.csv")
    print("the end")
    
    
# run main function 
main()


### end of get_wind_vectors_v2.py ###