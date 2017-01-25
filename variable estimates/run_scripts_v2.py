
################# run_scripts_v2.py ##################
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
# statistics_for_overlapping_zones.py script (provided by NOAA) is required for the batchRasterBufferIntersect function
# buffer_variables.py contains many custom functions called by runScripts.py
# constant_values.py conatins all modifiable input values (e.g. input files, folder locations)

############## import required modules ###############
import os
import buffer_variables
import multiprocessing
import arcpy
import constant_values as values
import gc
arcpy.env.overwriteOutput = True
import shutil
import time
import constant_values
############## end of module import ##################



########### helper functions ################

# create all of the buffer zones for all of the air monitor partitions 
# INPUTS:
#    airMonitorPartitions (string array) - list of full filepaths to air monitor partition shapefiles
def makeBufferZones(airMonitorPartitions):
    pool = multiprocessing.Pool()
    pool.map(BufferVariables.makeMultipleBuffers, airMonitorPartitions) # run make buffer zones on parallel processes
    pool.close()
    pool.join()    
    print ("completed making buffer zones for all partitions")
    del pool
### end of makeBufferZones ###

# get the filepaths for raster files to process
# INPUTS:
#    airMonitor (string) - filepath for the air monitor shapefile to create buffers for
#    rasterValues (string array) - array of raster filepaths to process
def determineRasterList(airMonitor, rasterValues):
    for fileName in values.RASTER_LIST:
        rasterValues.append(fileName)
    if(len(values.MOSAIC_RASTER_LIST) > 0):
        for variable in values.MOSAIC_RASTER_LIST:
            zone = BufferVariables.determineAirMonitorZone(airMonitor)
            mosaicFilename = BufferVariables.determineMosaicFile(variable,zone,values.RASTER_TYPE)
            rasterValues.append(mosaicFilename)        
            
# get the filepaths for polyline files to process
# INPUTS:
#   airMonitor (string) - filepath for the air monitor shapefile to create line buffers for
#   polyLineValues (string array) - array of polyline filepaths to process
def determinePolylineList(airMonitor, polyLineValues):
    for fileName in values.POLYLINE_LIST:
        polyLineValues.append(fileName)
    if(len(values.POLLYLINE_MOSAIC_LIST) > 0):
        for variable in values.POLLYLINE_MOSAIC_LIST:
            zone = BufferVariables.determineAirMonitorZone(airMonitor)
            mosaicFilename = BufferVariables.determineMosaicFile(variable,zone,values.POLYLINE_TYPE)
            polyLineValues.append(mosaicFilename) 
   
# get the filepaths for point files to process using buffers
# INPUTS:
#   airMonitor (string) - filepath for the air monitor shapefile to create point buffers for
#   pointBufferList (string array) - array of point filepaths to process
def determinePointBufferList(airMonitor, pointBufferList):
    for fileName in values.POINT_BUFFER_LIST:
        pointBufferList.append(fileName)
    

# get the filepaths for point files to process without using buffers
# INPUTS:
#    airMonitor (string) - filepath for the air monitor shapefile to extract point estimates for
#    pointList (string array) - array of point filepaths to process
def determinePointList(airMonitor, pointList):
    for file in values.POINT_LIST:
        pointList.append(file)
    if(len(values.POINT_MOSAIC_LIST) > 0):
        for variable in values.POINT_MOSAIC_LIST:
            zone = BufferVariables.determineAirMonitorZone(airMonitor)
            mosaicFilename = BufferVariables.determineMosaicFile(variable,zone,values.POINT_TYPE)
            pointList.append(mosaicFilename)

# setup and calculate average values for a buffer zone
# INPUTS:
#   partitionFolderOut (string) - filepath where results of buffer analyses will be stored in shapefile format
#   masterBufferFile (string array) - array of variable filepaths to extract buffer values from
#   identifier (int) - integer indicating which region of the world the air monitor partition is located in
#   bufferVal (int) - buffer distance in meters
#   airMonitor (string) - filepath to input air monitor shapefile containing centroid locations to create buffers around
#   fileList (list of string arrays) - custom object containing lists of both polyline and raster variables to create buffer variables for
def processBufferVariables(partitionFolderOut, masterBufferFile, identifier, bufferVal, airMonitor,variableType, fileList):
    readyToJoin = False
    continueVar = True
    while (continueVar):
        try:
            argumentList = []
            argumentList2 = []
            argumentList3 = []
            argVars = []
            if(variableType == values.PARALLEL_PROCESSING):
                rasterList = fileList[0]
                polylineList = fileList[1]  
                if(len(fileList)>0):
                    print("the rasterListLength is " + str(len(rasterList)))
                    pool = multiprocessing.Pool(len(rasterList))
                    argumentList = BufferVariables.createArgumentList(rasterList, partitionFolderOut, masterBufferFile,identifier,bufferVal,airMonitor)  
                    result = pool.map_async(BufferVariables.multi_run_raster_wrapper,argumentList) # calculate average raster values on parellel processors
                    argumentList2 = []
                    argumentList2 = BufferVariables.createArgumentList(polylineList, partitionFolderOut, masterBufferFile,identifier,bufferVal,airMonitor)  
                    pool2 = multiprocessing.Pool(len(polylineList))
                    result2 = pool2.map(BufferVariables.multi_run_polyline_wrapper,argumentList2)  # calculate average polyline values on parallel processors
                    pool2.close()
                    pool2.join()
                    BufferVariables.testProgress(result)
                    readyToJoin = BufferVariables.testFileCompletion(argumentList)  
            elif(variableType == values.RASTER_TYPE): # if the variable files are from rasters, run the raster wrapper function
                if(len(fileList)>0):
                    pool = multiprocessing.Pool(len(fileList))
                    argumentList = BufferVariables.createArgumentList(fileList, partitionFolderOut, masterBufferFile,identifier,bufferVal,airMonitor)  
                    result = pool.map_async(BufferVariables.multi_run_raster_wrapper,argumentList) # calculate average raster values on parellel processors
                    BufferVariables.testProgress(result)
                    readyToJoin = BufferVariables.testFileCompletion(argumentList)
            elif(variableType == values.POLYLINE_TYPE): # if the variable files are from polyline shp files, ru nthe polyline wrapper function      
                if(len(fileList)>0):
                    pool = multiprocessing.Pool(len(fileList))
                    argumentList2 = BufferVariables.createArgumentList(fileList, partitionFolderOut, masterBufferFile,identifier,bufferVal,airMonitor)  
                    result = pool.map(BufferVariables.multi_run_polyline_wrapper,argumentList2)  # calculate average polyline values on parallel processors
                    pool.close()
                    pool.join() 
                    readyToJoin=True
            elif(variableType==values.POINT_BUFFER_TYPE):
                if(len(fileList)>0):
                    argumentList3 = BufferVariables.createArgumentList(fileList, partitionFolderOut, masterBufferFile,identifier,bufferVal,airMonitor)
                    BufferVariables.pointBufferIntersect(argumentList3[0][0], argumentList3[0][1],argumentList3[0][2], argumentList3[0][3], argumentList3[0][4])
                    readyToJoin=True
            if(readyToJoin):
                for argument in argumentList: # for each variable that was used to calculate an average values, add the value to the air monitor partition shp file
                    BufferVariables.addVariableToPartition(argument, airMonitor,values.RASTER_TYPE)
                for argument in argumentList2:
                    BufferVariables.addVariableToPartition(argument, airMonitor,values.POLYLINE_TYPE)
                for argument in argumentList3:
                    BufferVariables.addVariableToPartition(argument, airMonitor,values.POINT_BUFFER_TYPE)
            else:
                raise Exception
            #print("sucessfully added buffer variables")
        except Exception as e:
            print("couldn't process variables, loop will cycle again" + str(e))
            pool.terminate()
            continue
        finally:
            if(variableType == values.RASTER_TYPE) or(variableType == values.PARALLEL_PROCESSING):
                try:
                    pool.terminate()
                except:
                    print("could not terminate pool")
                try:
                    dirList = os.listdir(values.RESULTS_FOLDER + values.TEMP_STATS_WORKSPACE)
                    for dirFolder in dirList:
                        path = values.RESULTS_FOLDER + values.TEMP_STATS_WORKSPACE + "/" + dirFolder
                        dirList2 = os.listdir(path)
                        for dirFolder2 in dirList2:
                            path2 = path + "/" + dirFolder2
                            dirList3 = os.listdir(path2)
                            for dirFolder3 in dirList3:
                                path3 = path2 + "/" + dirFolder3
                                try:
                                    shutil.rmtree(path3)    
                                except:
                                    exceptTemp = 1
                                try:
                                    fileList = os.listdir(path3)
                                    print(fileList)
                                    for filename in fileList:
                                        try:
                                            os.remove(path2 +"/" + filename)
                                        except:
                                            exceptTemp = 1
                                except:
                                    excpetTemp = 1               
                except:                            
                    exceptTemp = 1
        break
    print("finished calculating " + str(variableType) + " type values for buffer " + str(buffer) + " of air monitor partition " + identifier)     
    try:
        del pool, argVars, argumentList, fileList
    except:
        print("couldn't delete pool")
    
### end of processBufferVariables ###       
  
    
########## end of helper functions ###############
 
 
 
 
 
############ main function ##################
def main():
    print("running main")
    zonesDefined = BufferVariables.assignZones()
    if not os.path.exists(constant_values.RESULTS_FOLDER + constant_values.TEMP_STATS_WORKSPACE): os.makedirs(constant_values.RESULTS_FOLDER + constant_values.TEMP_STATS_WORKSPACE)   
    print("defined buffer zones")
    airMonitorPartitions = BufferVariables.partitionShapefile(zonesDefined) # partition air monitor stations
    print("defined air monitor partitions")
    airMonitorPartitions = airMonitorPartitions[7:len(airMonitorPartitions)]
    makeBufferZones(airMonitorPartitions) # make buffer zones for each partition
    i=0
    for airMonitor in airMonitorPartitions: # for each air monitor partition
        startTime = time.time()
        if(i>=0):
            rasterList = []
            polyLineList = []
            pointList = []
            pointBufferList = []
            determineRasterList(airMonitor,rasterList)
            determinePolylineList(airMonitor, polyLineList)
            determinePointList(airMonitor, pointList)
            determinePointBufferList(airMonitor, pointBufferList)
            print(pointList)
            #BufferVariables.runPointAnalysis(airMonitor, pointList)
            identifier = BufferVariables.determineAirMonitorIdentifier(airMonitor) # determine the partition number
            partitionFolderOut = values.RESULTS_FOLDER + values.KEYWORD + identifier + "/" 
            for buffer in values.BUFFER_DISTANCE: # for each buffer radius          
                bufferFilename = "buffer" + str(buffer) + "m.shp" 
                masterBufferFile = values.RESULTS_FOLDER + values.KEYWORD + identifier + values.BUFFER_EXTENSION + bufferFilename 
                maxThreads = multiprocessing.cpu_count()*2 -2
                if(maxThreads >= len(rasterList) + len(polyLineList)):
                    print("running polyline and raster buffer variables in parallel")
                    parallelList = [rasterList,polyLineList]
                    processBufferVariables(partitionFolderOut, masterBufferFile, identifier, buffer, airMonitor, values.PARALLEL_PROCESSING, parallelList)
                    processBufferVariables(partitionFolderOut, masterBufferFile, identifier, buffer, airMonitor, values.POINT_BUFFER_TYPE, pointBufferList)
                else:
                    processBufferVariables(partitionFolderOut, masterBufferFile, identifier, buffer, airMonitor,values.RASTER_TYPE, rasterList) # get average raster values
                    processBufferVariables(partitionFolderOut, masterBufferFile, identifier, buffer, airMonitor,values.POLYLINE_TYPE,polyLineList) # get average polyline values
                    processBufferVariables(partitionFolderOut, masterBufferFile, identifier, buffer, airMonitor,values.POINT_BUFFER_TYPE,pointBufferList) # get average point values
                print ("completed buffer distance " + str(buffer) + " for air Monitor partition " + str(identifier))
            print("time required to process partition " + str(identifier) + ": " + str(time.time()-startTime))
        i+=1
        print("completed gathering buffer values for all air monitoring station partitions")
    arcpy.Merge_management(inputs=airMonitorPartitions,output=values.RESULTS_FOLDER + "final.shp",field_mappings="#")  
    print ("completed running the main script")

### end of main function ###
    

# run the main function
if __name__ == '__main__':
    main()
    
    
########## end of run_scripts_v2.py ###############