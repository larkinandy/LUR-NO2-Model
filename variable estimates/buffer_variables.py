################# buffer_variables.py ##################
#
# Contains functions for automating several ArcGIS functions during the development of a Land Use Regression model.
# Functions include generating a series of buffers around points in a shapefile, determining polyline length within each unique buffer,
# and determining arverage raster values within the buffer.
#
# Author: Andrew Larkin
# Created for: Perry Hystad, Oregon State University
# Last modified: 01/25/2017
#
# Requirements:
# ArcGIS with spatial extension (tested on ArcGIS v. 10.3.1)
# ArcGIS comptabile version of python (tested with Python v. 2.7)
# Python integrated development environment is highly recommended but not required
# StatisticsForOverlappingZones.py script (provided by NOAA) is required for the batchRasterBufferIntersect function
# constantvalues.py conatins all modifiable input values (e.g. input files, folder locations)


############## import required modules ###############
import arcpy
from arcpy import sa
arcpy.CheckOutExtension("spatial")
arcpy.env.overwriteOutput = True
from arcpy import env
import statistics_for_overlapping_zones as overlap
import constant_values as values
import gc
import os
import math
import time
import sys
import shutil
import constant_values
############## end of module import ##################



################# functions ##################################

# determine the unique number assigned to the partition dataset in process
# INPUTS:
#   airMonitorPartitionFile (string) - filepath to the air monitor partition of interest
# OUTPUTS:
#   identifier (string) - component of the filepath that uniquely identifies the air monitor of interest
def determineAirMonitorIdentifier(airMonitorPartitionFile):
    startLocation = airMonitorPartitionFile.rfind(values.KEYWORD) + len(values.KEYWORD) 
    endLocation = airMonitorPartitionFile.find(".shp")
    identifier = airMonitorPartitionFile[startLocation:endLocation]
    print("this air monitor's identity is " + str(identifier))
    return identifier
### end of determineAirMonitorIdentifier

# determine which zone each air monitor is located in
def assignZones():
    #relpath = os.path.dirname(os.path.realpath(__file__))
    zoneDefined = values.RESULTS_FOLDER + "w_zones.shp"
    print("defined zoneDefined")
    arcpy.SpatialJoin_analysis(values.INPUT_FOLDER + values.MONITOR_FILE,values.INPUT_FOLDER + values.ZONE_DEFINITIONS,zoneDefined,"JOIN_ONE_TO_ONE","KEEP_ALL","#","CLOSEST","#","#")
    print("completed zone assignments")
    return zoneDefined


# calculate point values for air monitoring stations
# INPUTS:
#   airMonitorFile (string) - filepath to the air monitor shapefile of interest
#   pointList (string array) - names of the variables to extract point values from
def runPointAnalysis(airMonitorFile, pointList):
    #calculate distance from coast to air monitor location
    arcpy.Near_analysis(values.INPUT_FOLDER + values.MONITOR_FILE, values.INPUT_FOLDER + values.COAST_BOUNDARY)  
    rasterArgumentsList = ""
    # for each raster file used as point value data source, add it to the extract multivalues list and run extractMutliValues with arcpy
    if(len(pointList) >0):
        for pointFile in pointList:
            rasterArgument = values.INPUT_FOLDER + pointFile
            outputFile = values.RESULTS_FOLDER + "pointAnalysisTemp.shp"
            arcpy.sa.ExtractValuesToPoints(airMonitorFile, rasterArgument,values.RESULTS_FOLDER + "pointAnalysisTemp.shp", "NONE")
            arcpy.AddField_management(outputFile, pointFile[0:2], "DOUBLE", "", "")
            expression2 = "!" + "RASTERVALU" + "! * 1"
            arcpy.CalculateField_management(outputFile, pointFile[0:2], expression2, "PYTHON")
            arcpy.DeleteField_management(outputFile, "RASTERVALU")
            arcpy.CopyFeatures_management(outputFile, airMonitorFile)  
    else:
        print("there are no point variables to process")
    print("completed calculating point values for the air monitor data set")
### end runPointAnalysis ###

# test whether threads have completed processing buffers.  If threads have not made progress within alloted time defined in the 
# constant values file, kill the threads and try again
# INPUTS:
#    result (python object) - custom python object with results of whether a batch of threads have completed or not
def testProgress(result):
    print("we're testing progress")
    timeDif = 0
    while(timeDif <values.RATER_PROCESS_WAIT_TIME):
        try:
            result.get(timeout=values.RATER_PROCESS_WAIT_TIME)
            if(result.successful()==True):
                break            
        except Exception as e:
            print("timeout" + str(e))
        timeDif = time.time() - os.path.getmtime(values.RESULTS_FOLDER + values.TEMP_STATS_WORKSPACE + "/" + values.TEST_PROGRESS_FILE)
        print(timeDif)
    print("multiprocessing has either completed or stopped progress")
### end of testProgress ###

# determine which zone an air monitor is in based on the unique identifier included in the filename
# INPUTS:
#    airMonitorPartitionFile (string) - air monitor filepath that includes unique identifier
# OUTPUTS:
#    zoneIdentifier (int) - integer indicating which region of the world the air monitor is located in
def determineAirMonitorZone(airMonitorPartitionFile):
    startLocation = airMonitorPartitionFile.rfind(values.ZONE_KEYWORD)
    endLocation = airMonitorPartitionFile.rfind(values.PARTITION_KEYWORD)
    zoneIdentifier = airMonitorPartitionFile[startLocation+1:endLocation]
    return zoneIdentifier
### end of determineAirMonitorZone

# determine the exact raster from a mosaic should be used in a buffer analysis
# INPUTS:
#   mosaicFolder (string) - folder that contains all mosaic files for the variable of interest
#   zoneIdentifier (int) - integer identifying the region of world the air monitor is located in
#   fileType (string) - indicates what type of data (line, raster, etc) the mosaic of interest is
# OUTPUTS:
#   selectFile (string) - filepath for specific file within the mosaic to use for buffer analysis
def determineMosaicFile(mosaicFolder, zoneIdentifier, fileType):
    selectFile=[]
    fileList = os.listdir(values.INPUT_FOLDER + mosaicFolder)
    if fileType == values.RASTER_TYPE or fileType == values.POINT_TYPE:
        extension = ".tif"
    elif(fileType == values.POLYLINE_TYPE):
        extension = ".shp"
    fileList = os.listdir(values.INPUT_FOLDER + mosaicFolder)
    for file in fileList:
        startLocation = file.rfind("z")
        endLocation = file.rfind(extension)
        zoneValue = file[startLocation+1:endLocation]
        if zoneValue == zoneIdentifier:
            selectFile = mosaicFolder + "/" + file[0:endLocation] + extension 
    return selectFile

# partition a shapefile with a large number of air monitoring stations into several shapefiles with a smaller
# number of monitoring stations per file is determined by the PARTITION_SIZE variable in the constant_values.py file
# INPUTS:
#   airMonitorFile (string) - filepath to the shapefile containing all air monitor locations
# OUTPUTS:
#   shapeFileList (string array) - filepaths to all of the air monitor partitions.
def partitionShapefile(airMonitorFile):
    newPath = values.RESULTS_FOLDER
    shapeFileList = []
    with arcpy.da.SearchCursor(airMonitorFile,"zone") as cursor:      
        for zone in sorted({row[0] for row in cursor}):
            whereClause = '"zone" = ' + str(zone) 
            print(zone)
            tempFile = newPath + "spat_join_temp.shp"
            lyr = arcpy.MakeFeatureLayer_management(airMonitorFile)
            arcpy.Select_analysis(lyr, tempFile, whereClause)
            getCount = int(arcpy.GetCount_management(tempFile).getOutput(0))
            numPartitions = int(math.ceil(float(getCount)/values.PARTITION_SIZE))
            print("the number of partitions for zone " + str(zone) + "is " + str(numPartitions))
            for i in range(0,numPartitions):
                newPath = values.RESULTS_FOLDER + "/" + values.KEYWORD + "z" + str(zone) + "i" + str(i)
                if not os.path.exists(newPath): os.makedirs(newPath)                 
                outputFileName =newPath + "/" + values.KEYWORD +  "z" + str(zone) + "i" + str(i) + ".shp"
                whereClause = '"FID" >= ' + str(i*values.PARTITION_SIZE) + ' AND "FID" < ' + str(min(getCount+1, (i+1)*values.PARTITION_SIZE))
                print(whereClause)
                arcpy.Select_analysis(tempFile, outputFileName, whereClause)
                print("created partition zone " + str(zone) + ", id " +  str(i))
                shapeFileList.append(outputFileName)
            print("completed zone " + str(zone))
    print("finished paritionShapeFile")
    print(shapeFileList)    
    return(shapeFileList)




# add varaible values to the partition air monitor shapefile
# INPUTS:
#    argument (string) - name of the variable to add to the attribute table
#    airMonitor (string) - filepath to the airMonitor shapefile to add the attribute to
#    valueType (int) - type of variable that corresponds to added value (polyline buffer, raster buffer, point, etc)
def addVariableToPartition(argument, airMonitor,valueType):
    bufferIndex = str(argument[0]).rfind("buffer") + 6
    variableIdent = argument[0][bufferIndex:-4]
    if(valueType == values.RASTER_TYPE):
        valueFile = argument[2] + variableIdent + ".shp"
        arcpy.JoinField_management(airMonitor,"FID",valueFile,values.AIRMONITOR_ID,variableIdent)
    elif(valueType == values.POLYLINE_TYPE):
        valueFile = argument[2] + variableIdent + "d.shp"
        arcpy.JoinField_management(airMonitor,"FID",valueFile,values.BUFFER_ID,variableIdent)
    elif(valueType == values.POINT_BUFFER_TYPE):
        valueFile = argument[2] + variableIdent + "d.shp"
        arcpy.JoinField_management(airMonitor,"FID",valueFile,values.BUFFER_ID,variableIdent)        
    print ("completed adding the " + variableIdent + " to the air monitor partition " + airMonitor)
### end of addVariableToPartition


# make a buffer of given input distance from the points in the given input file.  
# write the results into the designated outputFolder
# INPUTS:
#    bufferDistance (int) - buffer radius, in meters
#    inputFile (string) - filepath to the shapefile with locations to create buffers for
#    dataFolderOut (string) - filepath to folder that will contain resulting buffer shapefile
# OUTPUTS:
#    outputFilepath (string) - filepath to the buffer shapefile that was created
def makeBuffer(bufferDistance, inputFile, dataFolderOut):  
    # for each bufferDistance to be created, run through the ArcGIS buffer analysis tool and create an output .shp file 
    inputFilePath = inputFile
    outputFile = "buffer" + str(bufferDistance) + "m.shp"
    outputFilePath = dataFolderOut + outputFile
    buffDistance = str(bufferDistance) + " Meters"
    if(os.path.exists(outputFilePath)):
        print("warning: " + outputFilePath + "already exists. File will be overwritten")
    arcpy.Buffer_analysis(in_features=inputFilePath,out_feature_class=outputFilePath,buffer_distance_or_field=buffDistance,
                          line_side="FULL",line_end_type="ROUND",dissolve_option="NONE",dissolve_field="#")
    return outputFilePath 
### end of makeBuffer ### 


# create buffer shapefiles and a folder to store them in
# INPUTS:
#    partitionFile (string) - filepath to the air monitor partition to create buffers for
def makeMultipleBuffers(partitionFile):
    startLocation = partitionFile.rfind(values.KEYWORD) + len(values.KEYWORD)
    endLocation = partitionFile.find(".shp")
    identifier = partitionFile[startLocation:endLocation]    
    partitionFolderOut = values.RESULTS_FOLDER  + values.KEYWORD + identifier
    bufferFolder = partitionFolderOut + values.BUFFER_EXTENSION
    if not os.path.exists(bufferFolder): os.makedirs(bufferFolder)
    for buffer in values.BUFFER_DISTANCE:
        makeBuffer(buffer,partitionFile, bufferFolder)
### end of makeBuffers

# create a copy of the buffer file for each parallel processing thread
# INPUTS:
#   masterBufferFile (string) - file containing template buffers to copy
#   variable (string array) - array of unique identifiers for each thread
#   bufferVal (string) - buffer radius, in meters
#   partitionFolderOut (string) - filepath where the copy of the buffer will be created
# OUPUTS:
#   tempBufferFile (string) - filepath to the copy of the buffer file
def createBufferFileCopy(masterBufferFile, variable, identifier, bufferVal, partitionFolderOut):
    tempFileName = "buffer" + variable[0] + variable[1] + bufferVal + "m.shp"
    tempBufferFile = partitionFolderOut + values.BUFFER_EXTENSION + tempFileName
    arcpy.CopyFeatures_management(masterBufferFile, tempBufferFile)
    return tempBufferFile
### end of createBufferFileCopy


# create a list of arguments to pass into the wrapper function for parallel processing of a variable
# INPUTS:
#   variableList (string array) - list of variables to create buffer estimates for
#   partitionFolderout (string) - folder where buffer files for each variable will be created
#   masterBufferFile (string) - filepath to the original buffer file that will be copied for each variable
#   identifier (string) - number that identifies which region of the world the buffer files are located in
#   bufferVal (int) - buffer radius, in meters
#   airMonitor (string) - filepath to the air monitor shapefile location
# OUTPUTS:
#   argumentList (obj array) - array of custom object containing arguments needed to spawn mupltiple threads, one thread
#                              for each variable
def createArgumentList(variableList, partitionFolderOut, masterBufferFile,identifier,bufferVal,airMonitor):
    argumentList = []
    for variable in variableList: # for each variable that we want to proces (i.e. included in the list)
        variableIdent = variable[0] + variable[1]
        variableOutputFolder = partitionFolderOut + variable[0] + variable[1] + "/" 
        tempBufferFile = createBufferFileCopy(masterBufferFile, variable, identifier, str(bufferVal), partitionFolderOut)
        # add the created list of arguments to the master list of arguments to be used in the call for parallel processing
        argumentList.append([tempBufferFile, values.INPUT_FOLDER + variable, variableOutputFolder, airMonitor, str(bufferVal)]) 
    print("completed making arguments for parallel processing of buffer varaible analysis")
    return argumentList
### end of createArguemntList


# calculate the length of a polyline file within given buffer areas
# INPUTS:
#   bufferFile (string) - filepath to the shapefle that describes the buffer areas of interest
#   variableFile (string) - filepath to the shapefile that contains the polyline data of interest
#   varibleOutputFolder (string) - folderpath to where temporary shapefiles are stored
#   airMonitorFile (string) - filepath to air monitor locations
#   bufferSize (int) - radius of the buffer, in meters
def polylineBufferIntersect(bufferFile, variableFile,  variableOutputFolder, airMonitorFile, bufferSize):
    bufferIndex = str(bufferFile).rfind("buffer") + 6
    variableIndex = str(variableFile).rfind("/") + 1
    variableIdent = bufferFile[bufferIndex:-4]
    intersectFile = variableOutputFolder + variableIdent + ".shp"
    dissolveFile = variableOutputFolder + variableIdent + "d.shp"
    if not os.path.exists(variableOutputFolder): os.makedirs(variableOutputFolder)   
    tableFile = variableOutputFolder + variableIdent + ".shp"
    if(os.path.isfile(tableFile)):
        return    
    fieldName = variableIdent # new field to add to feature class table 
    
    # intersect the buffer
    arcpy.Intersect_analysis(bufferFile +";" + variableFile, intersectFile, "ALL", "", "INPUT")
   
    # combine line segments within each buffer zone (overlapping areas are accounted for)
    arcpy.Dissolve_management(intersectFile,dissolveFile,values.BUFFER_ID,"#","MULTI_PART","DISSOLVE_LINES") 
   
    # create a new field in the dissolved feature class
    arcpy.AddField_management(dissolveFile, fieldName, "DOUBLE", "", "")
   
    # calculate the total length of polyline segments in each unique buffer and store the result in the new field
    arcpy.CalculateField_management(dissolveFile, fieldName, values.LENGTH_COMMAND, "PYTHON") 
    
    del (bufferFile, variableFile, airMonitorFile, bufferSize, bufferIndex, variableIndex, 
         variableIdent, variableOutputFolder, intersectFile, dissolveFile, tableFile, fieldName)
    gc.collect()    

# perform buffer analysis on point variables
# INPUTS:
#   bufferFile (string) - filepath to the template buffer shapefile
#   variableFile (string) - folder to the buffer shapefile copy created for the analysis
#   airMonitorFile (string) - filepath to the shapefile describing air monitor locations
#   bufferSize (int) - radius of the buffer, in meters
def pointBufferIntersect(bufferFile, variableFile,  variableOutputFolder, airMonitorFile, bufferSize):
    bufferIndex = str(bufferFile).rfind("buffer") + 6
    variableIndex = str(variableFile).rfind("/") + 1
    variableIdent = bufferFile[bufferIndex:-4]
    intersectFile = variableOutputFolder + variableIdent + ".shp"
    dissolveFile = variableOutputFolder + variableIdent + "d.shp"    
    if not os.path.exists(variableOutputFolder): os.makedirs(variableOutputFolder)   
    tableFile = variableOutputFolder + variableIdent + ".shp"
    if(os.path.isfile(tableFile)):
        return    
    fieldName = variableIdent # new field to add to feature class table 
    
    # intersect the buffer
    arcpy.Intersect_analysis(bufferFile +";" + variableFile, intersectFile, "ALL", "", "INPUT")
    
    # combine line segments within each buffer zone (overlapping areas are accounted for)
    # Replace a layer/table view name with a path to a dataset (which can be a layer file) or create the layer/table view within the script
    # The following inputs are layers or table views: "powerPlants_Intersect"
    arcpy.Dissolve_management(intersectFile,dissolveFile,values.BUFFER_ID,"carbon_200 SUM","MULTI_PART","DISSOLVE_LINES")    
     
    # create a new field in the dissolved feature class
    arcpy.AddField_management(dissolveFile, fieldName, "DOUBLE", "", "")
    
    # calculate the total length of polyline segments in each unique buffer and store the result in the new field
    arcpy.CalculateField_management(dissolveFile, fieldName, values.CARBON_COMMAND, "PYTHON")   
    
    del (bufferFile, variableFile, airMonitorFile, bufferIndex, variableIndex, 
         variableIdent, variableOutputFolder, intersectFile, dissolveFile, tableFile, fieldName)
    gc.collect()    

# determine whether all buffers have been made
# INPUTS:
#   arugmentList (array obj) - array of unique python objects, used to spawn threads for multivariate analysis
# OUTPUTS:
#   fileCompletion (boolean) - whether or not all buffer distances have completed 
def testFileCompletion(argumentList):
    fileCompletion = True
    for argument in argumentList:
        bufferIndex = str(argument[0]).rfind("buffer") + 6
        variableIdent = argument[0][bufferIndex:-4]
        valueFile = argument[2] + variableIdent + ".shp"
        testTemp = 0
        if(os.path.isfile(valueFile)):
            for field in arcpy.ListFields(valueFile):
                if(field.name == argument[2] + variableIdent):
                    testTemp = 1
                else:
                    testTemp = max(testTemp,0)
        fileCompletion = fileCompletion * os.path.isfile(valueFile)
    return(fileCompletion)

# parallel processing natively only accepts a single argument for input.  
# this is a wrapper function to accept a single list of inputs from the parall
def multi_run_raster_wrapper(args):
    return rasterBufferIntersect(*args)
### end of multi_run_wrapper


# parallel processing natively only accepts a single argument for input.  
# this is a wrapper function to accept a single list of inputs from the parall
def multi_run_polyline_wrapper(args):
    return polylineBufferIntersect(*args)
### end of multi_run_wrapper

# calculate average values of a raster within various buffer zones in a shapefile
# INPUTS:
#   bufferFile (string) - filepath to template buffer shapefile
#   varaibleFile (string) - filepath to the variable input dataset
#   varaibleOutputFolder (string) - folderpath to where intermediate results will be stored
#   airMonitorFile (string) - filepath to the shapefile describing air monitor locations
#   bufferSize (int) - radius of buffer, in meters
def rasterBufferIntersect(bufferFile, variableFile, variableOutputFolder, airMonitorFile, bufferSize):
    continueVar = 1
    var1 = variableFile
    bufferIndex = str(bufferFile).rfind("buffer") + 6
    variableIndex = str(variableFile).rfind("/") + 1
    variableIdent = bufferFile[bufferIndex:-4]
    partitionIdentStart = variableOutputFolder.rfind("Partition") + len("Partition")
    partitionIdentEnd = variableOutputFolder.rfind("/")
    partitionId = variableOutputFolder[partitionIdentStart:partitionIdentEnd]
    if not os.path.exists(variableOutputFolder): os.makedirs(variableOutputFolder)
    tableFile = variableOutputFolder + variableIdent + ".shp"   
    tempTestValue = 1
    while(continueVar==1):
        if(os.path.isfile(tableFile)):
            print(str(tableFile) + " already exists")
            sys.stdout.flush()
            for field in arcpy.ListFields(tableFile):
                sys.stdout.flush()
                if field.name == variableIdent:
                    continueVar = 0
                    tempTestValue = 0
                    print("the field already exists, no need to rerun for variable " + str(variableIdent))
                    sys.stdout.flush()
                    return
            if(tempTestValue == 1):
                try:
                    os.remove(variableOutputFolder + variableIdent + ".shp")
                    os.remove(variableOutputFolder + variableIdent + ".cpg")
                    os.remove(variableOutputFolder + variableIdent + ".shx")
                    os.remove(variableOutputFolder + variableIdent + ".dbf")
                    os.remove(variableOutputFolder + variableIdent + ".prj")
                    #deleteTest = arcpy.Delete_management(tableFile)
                    print("deleted corrupt shapefile " + tableFile )
                except exception as e:
                    print("cannot delete shapefile " + str(tableFile) + str(e))
                sys.stdout.flush()
        if(continueVar==1):
            try:     
                raster = arcpy.Raster(var1)
                arcpy.env.cellSize = raster.meanCellHeight/10.0
                expression2 = []
                try:
                    continueVar = overlap.runFeatureRasterIntersect(bufferFile, values.TABLE_ID, variableFile, tableFile, variableFile[variableIndex:variableIndex+2],partitionId,str(bufferSize))    
                except:
                    print("could not run overlap.runfeatureIntersect")
                try:
                    for field in arcpy.ListFields(tableFile, "*_MEAN"):
                        meanField = field.name
                        expression2 = "!" + meanField + "! * 1"
                except:
                    expression2 = "-9999"
                    print("-9999 was used")
                arcpy.AddField_management(tableFile, variableIdent, "DOUBLE", "", "")  
                arcpy.CalculateField_management(tableFile, variableIdent, expression2, "PYTHON")                 
                print(variableIdent + " raster for buffer size " + str(bufferSize) + " completed")
                del (bufferFile, variableFile, tableFile, raster, airMonitorFile, bufferSize)
                gc.collect()
                print("finished raster analysis")
                try:
                    tempWorkspace = constant_values.RESULTS_FOLDER + constant_values.TEMP_STATS_WORKSPACE + "/" + variableFile[variableIndex:variableIndex+2] + str(partitionId) + str(bufferSize) + "/zonalStats"
                    shutil.rmtree(tempWorkspace)
                except:
                    print("could not remove the tempworkspace")                
            except Exception as e:
                print("raster analysis failed " + str(e) + " " + variableIdent)
                continueVar = 1
                continue
        try:
            sys.stdout.flush()
        except:
            print("flush didn't work")
        continueVar = 0
### end of rasterBufferIntersect ###  
  
################# end of functions ############################


############### end of buffer_variables.py ###############