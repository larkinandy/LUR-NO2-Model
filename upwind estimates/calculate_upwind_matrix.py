#############################################
# calculate_upwind_matrix.py 
# Author: Andrew Larkin
# Developed for Perry Hystad, Oregon State University
# Date last modified: Jan 25, 2017
# Description: This script combines wind direction measurements from a csv file with
# air monitor shapefiles and buffered major and minor roads around air monitors
# to estimate the amount of roads upwind and downwind from the air monitor in a given year

# import libraries
import csv
import math
import arcpy
import os
import sys

arcpy.env.overwriteOutput=True      #  overwrite to overwrite temp.shp and temp2.shp


parentFolder = "insert parent folder here" # main folder containing all files of interest
testInputFolder = parentFolder + "calculateWindVectors/dallasExample/"    # folder where all csv files and shapefiles are located
inputShapeFolder = testInputFolder + "inputShapefiles/" # folder where shapefiles are located
windFilename = "shpfl.shp"                              # name for the temporary wind vector shapefile
inputCSVFilepath = testInputFolder + "dallasWind.csv"           # CSV containing wind direction and date stamps
CSV_DICT = ["latitude","longitude"]     # values are stored in latitude and longitude
WIND_FIELD = "UpWind"                      # attribute field in temp and temp2 files
LATITUDE_FIELD = "latitude"                 # name of the attribute field in the shapefile containing latitude coords
LONGITUDE_FIELD = "longitude"               # name of the attirbute field in the shapefile containing longitude coords
NUM_TIME_SAMPLES = 365*4

sr = arcpy.SpatialReference()   # WGS84 spatial reference
sr.factoryCode = 4326
sr.create()

# the following components are used for the Haversine formula.  For more details go to
# http://www.movable-type.co.uk/scripts/latlong.html
# https://www.eol.ucar.edu/content/wind-direction-quick-reference

EARTH_RADIUS = 6378140.0                # used to calculate wind direction
distance = 55000.0/EARTH_RADIUS         # used to calculate bearing and wind direction


################### helper functions #######################


# create a 360 element array, with one element for each angular degree
def createAngularMatrix():
    angMatrix = [0] *360
    return(angMatrix)



# determine how often each of the 360 angular degrees are upwind from an air monitor.  For this calculation,
# upwind is defined as being within the windDegree number of degrees from directly upwind
# INPUTS:
#    windDirection (int) - direction of the wind for a single 6 hour sample
#    angMatrix (int array) - array that keeps track of the number of 6 hour samples are upwind for each 1 degree division of angle
def addWindExposure(windDirection,angMatrix):
    windDegree = 23     # tolerance for how many degrees deviation is allowed from directly upwind to still be classified as upwind
    
    if(windDirection == 360):        # if the wind value is 360, change to 0, which is the angular equivalent of 360
        windDirection = 0
    
    if(windDirection <windDegree):   # if the wind direction is less than the windDegree value, some values less than 0 (i.e. around 350 or so) are also upwind
        runWindComp(0,windDirection+windDegree,angMatrix) # calculate the number of days the degree angles to the right of due north are upwind
        leftSideMin = 360+windDirection-windDegree  # determine the minimum angle to the left of due north that is upwind
        runWindComp(leftSideMin,360,angMatrix) # calculate number of days the degree angles to the left of due north are upwind

    elif(windDirection > 338):       # if the wind direction is slightly to the left of due north, some values to the right of due north are also upwind
        runWindComp(windDirection-windDegree,360,angMatrix) # calculate the number of days the degree values to the left of due north are upwind
        rightSideMax = 360-windDirection + windDegree # determine the maximum angle to the right of due north that is upwind
        runWindComp(0, rightSideMax,angMatrix) # calculate the number of days the degree angles to the right of due north are upwind

    else:
        runWindComp(windDirection-windDegree,windDirection+windDegree,angMatrix) #calculate the number of days degree angles are upwind
            
            
# for each angle element in angMatrix, increment if they are upwind for the given wind direction            
# INPUTS:
#   windMinDegree (int) - minimum angle that is within error to be classified as upwind for a single measurement
#   windMaxDegree (int) - maximum angle that is within error to be classified as upwind for a single measurement
#   angMatrix (int array) - number of 6 hour samples that each angular degree from 0 to 359 is classified as upwind
def runWindComp(windMinDegree,windMaxDegree,angMatrix):
    
    for angle in range(0,360):
        if(angle >= windMinDegree and angle <= windMaxDegree):
            angMatrix[angle]+=1


# find the index in the two input datasets where the latitude and longitude values match
# INPUTS:
#   infoArr (obj arr) - custom data structure containing data for each row in a shapefile
#   latLongArr (obj arr) - custom data structure containing lat/long coordinates from a csv file
def matchLatLong(infoArr,latLongArr):
    i=0
    for latVal in infoArr[0]:
        if(latVal == latLongArr[0] and infoArr[1][i] == latLongArr[1]):
            return(i)
        i+=1
            

# read in a csv file and create a shapefile with angular degrees containing number of days each angle is upwind
# INPUTS:
#   inputFile (string) - filepath to the csv file to load
#   infoArr (obj array) - custom object containing shapefile meta information such as lat/long for each row in the shapefile
# OUTPUTS:
#   angMatrix (int array) - array in which each element in the array is the number of 6 hour samples an angular degree is upwind from centroid
def readCSVFile(inputFile,infoArr):
    latLongArray = [0,0]
    i=0
    
    with open(inputFile) as f: # open inthe input file  
        reader = csv.DictReader(f, dialect=f) # identify the names of the columns (i.e. variables) in the csv file

        for row in reader: # read a row as {column1: value1, column2: value2,...}
            wndArr = []     # array to hold number of days each angle is upwind
           
            for (k,v) in row.items(): # go over each column name and value 
              
                if(k in CSV_DICT): # if the variable name is latitude or longitude
                    latLongArray[CSV_DICT.index(k)] = float(v)      # store the value in the latitude longitude paired arra
                else:
                    wndArr.append(float(v))                         # if a wind direction value, store in the wind array

            indexNum = matchLatLong(infoArr,latLongArray)           # find the index in the the list of shapefiles that matches the current row in the csv file
            filename = infoArr[3][indexNum] + windFilename          # create a filename for a temporary wind direction shapefile, and create the shapefile
            arcpy.CreateFeatureclass_management(infoArr[3][indexNum], windFilename,"POLYGON",'#','#','#',sr)
            arcpy.AddField_management(filename, WIND_FIELD, "Short")   
            j=0
            angMatrix = createAngularMatrix()       # create a 360 element array, one element for each angular degree
            for wndDrct in wndArr:
                addWindExposure(wndDrct,angMatrix)   # calculate the number of 6 hour samples in which each angular degree is upwind
            for angle in angMatrix:
                calcTriangle(j+0.5,latLongArray[0],latLongArray[1],distance,filename,angle) # add 360 polygons to wind shapefile, one polygon for each angular degree
                j+=1
            i+=1
        del reader, row   
        return(angMatrix)       # retrun array with number of time points each angular degree is upwind


# calculate latitude and longitude coordinates on the outside of the buffer region
# INPUTS:
#   origLat (float) - latitude coordinate
#   origLong (float) - longitude coordinate
#   bearing (float) - angular direction relative to true north
#   distance (float) - distance that new coordinates should be from the center
# OUTPUTS:
#   newLatit (float) - latitude coordinate of outer buffer edge
#   newLongit (float) - longitude coordinate of outer buffer edge
def calcCoords(origLat,origLong,bearing,distance):    
    latRadians = origLat*0.0174533          # convert latitude from degrees to radians
    longRadians = origLong*0.0174533        # convert longitude from degrees to radians
    
    # calculate the new latitude coordinate for a point that is 'distance' units away using the haversine formula (see weblink above)
    newLatit = math.asin(math.sin(latRadians)*math.cos(distance) + math.cos(latRadians)*math.sin(distance)*math.cos(bearing))/0.0174533
    
    # calculate the new longitude coordinate for a point that is 'distance' units away using the haversine formula
    longitP1 = math.sin(bearing)*math.sin(distance)*math.cos(latRadians)
    longitP2 = math.cos(distance) - math.sin(latRadians)*math.sin(newLatit*0.0174533)    
    newLongit = longRadians + math.atan2(longitP1,longitP2)
    newLongit = newLongit/0.0174533
    newLong = (newLongit+540)%360-180
    return([newLatit,newLongit])


# create a triangle polygon in arcgis using 3 coordinate points.  The triangle will represent one angular degree of coverage
# INPUTS:
#   wndDrt (float) - angular wind direction, with 0 degrees corresponding to true north
#   strtLat (float) - latitude of the centroid
#   strtLong (float) - longitude of the centroid
#   inputFilename (string) - filepath to the shapefile to which the triangle will be added
#   timestamp (int) - time that corresponds to the wind direction, in epoch format
def calcTriangle(wndDrt,strtLat, strtLong,distance,inputFilename,timestamp):
    bearing1 = ((wndDrt + 180.5)%360.0)*0.0174533       # center of triangle + 0.5 degrees
    bearing2 = ((wndDrt - 180.5)%360.0)*0.0174533       # center of triangle - 0.5 degrees
    airMonitorLoc = arcpy.Point(strtLong,strtLat)       # defnie the air monitor location as one point of the triangle
    coordP1 = calcCoords(strtLat,strtLong,bearing1,distance)  # create another coordinate point
    upPoint1 = arcpy.Point(coordP1[1],coordP1[0])
    coordP2 = calcCoords(strtLat,strtLong,bearing2,distance)  # create another coordinate point
    upPoint2 = arcpy.Point(coordP2[1],coordP2[0])
    
    # Create a polygon geometry
    array = arcpy.Array([airMonitorLoc,upPoint1,upPoint2])
    polygon = arcpy.Polygon(array)
    # Open an InsertCursor and insert the new geometry
    cursor = arcpy.da.InsertCursor(inputFilename, ['SHAPE@', WIND_FIELD])
    cursor.insertRow([polygon,timestamp])
    # Delete cursor object
    del cursor    

# create wind triangles for all angular degrees (360 in total)
# INPUTS:
#   csvArray (float matrix) - matrix of latitude, longitude coordinates and various wind direction
#   inputFilename (string) - filepath to the shapefile where the wind triangles will be added
def calcWindTraingles(csvArray,inputFilename):
    
    latitude = csvArray[CSV_DICT.index("latitude")][0]
    longitude = csvArray[CSV_DICT.index("longitude")][0]
    
    for wndDrct in csvArray[CSV_DICT.index("windDrct")]:
        calcTriangle(wndDrct,latitude,longitude,distance,inputFilename)


# get latitude and longitude values from a shapefile
# INPUTS:
#   inputFilepath (string) - filepath for the shapefile of interest
# OUTPUTS:
#   latVal (float) - latitude coordinate
#   longVal (float) - longitude coordinate
def getShpflLatLong(inputFilepath):
    # Create a cursor on a feature class
    cur = arcpy.UpdateCursor(inputFilepath)
    # Loop through the rows in the attribute table
    for row in cur:
        latVal = row.getValue(LATITUDE_FIELD)
        longVal = row.getValue(LONGITUDE_FIELD)
        cur.updateRow(row) 
    del row, cur
    return([latVal,longVal])


# get latitude, longitude, folder name and filepath for all air monitor shapefiles
# INPUTS:
#   inputShapeFolder (string) - filepath to the input shapefolder
# OUTPUTS:
#   infoArr (obj arr) - custom data structure with latitude, longitude, filepath, and folder path for the shapefile of interest
def getShapeMeta(inputShapeFolder):
    filesToProcess = os.listdir(inputShapeFolder) # get the list of shapefiles to process
    latArr = []             # array to store latitude shapefile coords
    longArr = []            # array to store longitude shapefile coords 
    filepathArr = []        # array to store shpaefile filepaths
    folderpathArr = []      # array to sotre paths to shapefile folders
    
    # for each file in inputShapeFolder, determine if the file has .shp extension.  If so, get the latitude and longitude values, 
    # and record the filepath and folder path
    for filename in filesToProcess:    
        filepath = inputShapeFolder + filename + "/" + filename + ".shp"
        [latVal,longVal] = getShpflLatLong(filepath)
        latArr.append(latVal)
        longArr.append(longVal)
        folderpathArr.append(inputShapeFolder + filename + "/" )
        filepathArr.append(filepath)

    infoArr = [latArr,longArr,filepathArr,folderpathArr]        # combine all variables into a single array
    return(infoArr) 



# calculate average amount of roads for a single buffer distance and a single road type (major or minor roads) for a single air monitor
# INPUTS:
#   buffFolder (string) - folderpath containing files with road buffer data)
#   buffFile (string) - filename of buffer shapefile to process
#   airMonitorFolder (string) - folderpath containing air monitor shapefile of interest
#   airMonitorShp (string) - filename of air monitor shapefile of interest
def calcBufExpSingleBuff(buffFolder, buffFile,airMonitorFolder,airMonitorShp):
    
    inputFiles = [buffFolder + buffFile,airMonitorFolder+windFilename]          # define input files for arcpy intersect analysis
    arcpy.Intersect_analysis(inputFiles,airMonitorFolder + "temp.shp",'#','#','LINE') # run intersect between angular degree shapefile and buffered roads
    statVal = [["UpWind",'MEAN']]                         # define statistics to collect during the dissolve function
    # combine all roads within each angular degree and claculate summary statistics
    arcpy.Dissolve_management(airMonitorFolder + "temp.shp",airMonitorFolder + "temp2.shp","FID_shpFl",statVal,'#','DISSOLVE_LINES')   
    arcpy.AddField_management(airMonitorFolder + "temp2.shp", "AVG_LEN", "DOUBLE", "", "")
    avgLengthCommand = "!shape.length@kilometers! /  " + str(NUM_TIME_SAMPLES)
    arcpy.CalculateField_management(airMonitorFolder + "temp2.shp", "AVG_LEN", avgLengthCommand, "PYTHON")      
    fields = ['AVG_LEN','MEAN_UpWin']                    # define attribute fields in the dissolved results shapefile that contain results of interest
    sumVal = sumVals(airMonitorFolder + "temp2.shp",fields) # calculate the sum of roads for each anuglar degree, weighted by the number of 6hour intervals upwind
    fieldName = "u" + buffFile[0:len(buffFile)-4]            # define the name of the attribute field that will be written to the final output shapefile      
    arcpy.AddField_management(airMonitorShp,fieldName,"DOUBLE") 
    setVal(airMonitorShp,fieldName,sumVal)    # add the attribute field to the final shapefile  

# calculate average amount of roads for all buffer distances and road types (major and minor roads) for all air monitors
# INPUTS:
#   infoArr (obj array) - custom data structure containing meta data related to all air monitors in the master air monitor shapefile
def calcBufExp(infoArr):
    i=0
    
    for airMonitorFolder in infoArr[3]:         # for each air monitor shapefile
        for filename in os.listdir(airMonitorFolder + "mj/"):       # for each buffered major road shapefile
            if(filename[len(filename)-3:len(filename)]=="shp" and filename[len(filename)-5:len(filename)]!="d.shp"):
                calcBufExpSingleBuff(airMonitorFolder + "mj/",filename,airMonitorFolder,infoArr[2][i])      # claculate average major roads upwind
     
        for filename in os.listdir(airMonitorFolder + "mi/"):       # for each buffered minor road shapefile
            if(filename[len(filename)-3:len(filename)]=="shp" and filename[len(filename)-5:len(filename)]!="d.shp"):
                calcBufExpSingleBuff(airMonitorFolder + "mi/",filename,airMonitorFolder,infoArr[2][i])    # calculate average minor roads upwind
        i+=1
                   
    print("completed calcBufferExposure")


# add the average amount of road upwind for a buffered distance to the designated shapefile
# INPUTS:
#   inputFile (string) - filepath to the angular shapefile that should be updated
#   fieldName (string) - name of the attribute in the attribute table that should be updated
#   inVal (int) - updated value
def setVal(inputFile,fieldName,inVal):
    cur = arcpy.UpdateCursor(inputFile)
    for row in cur:
        row.setValue(fieldName, inVal)
        cur.updateRow(row)
    del row, cur



# get the average length of roads upwind from air monitor
# INPUTS:
#   inputFile (string) - filepath for the input file of interest
#   fields (string) - name of the attribute fields containing values needed to calculate length
# OUTPUTS:
#   sumVal (float) - length of roads upwind from air monitor
def sumVals(inputFile,fields):
    sumVal = 0
    numVals = 0
    with arcpy.da.SearchCursor(inputFile, fields) as cursor:
        for row in cursor:
            sumVal += row[0]*row[1]
            numVals +=1
        if(numVals>0):
            del row
        del cursor
    return(sumVal)


###################### main function ####################


def main():
    infoArr = getShapeMeta(inputShapeFolder)
    readCSVFile(inputCSVFilepath,infoArr)
    calcBufExp(infoArr)
    print("completed the main function")
    
    
    
main()




### end of calculate_upwind_matrix.py ###