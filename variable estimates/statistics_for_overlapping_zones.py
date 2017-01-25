#  Calculate Zonal Statistics for Overlapping Zones
# ---------------------------------------------------------------------------
# Import system modules
# Script curteousy of NOAA

# modified by Andrew Larkin for Perry Hystad, Oregon State University
# date last modified: 01/25/2017

import sys, string, os, gc, shutil
import arcpy
from arcpy import env
from arcpy.sa import *
import shutil
import gc
import time
import datetime
import random
import constantValues
import sys




# calculate buffer variables for all air monitors within a shapefile
# INPUTS:
#   inputFile1 (string) - filepath to the buffer shapefile
#   zoneField (int) - integer indicating region of the world the buffers are located within
#   outputFile (string) - filepath to where results will be written to in shapefile format
#   variableIdentity (string) - indicating whether a variable is a raster or a polyline file
#   partitionId (int) - unique Id for the current shapefile partition
#   bufferSize (int) - size of the buffer, in meters
def runFeatureRasterIntersect(inputFile1, zoneField, valueRaster, outputFile, variableIdentity, partitionId,bufferSize):
  continueVar = 1
  try: 
    scratchWorkspace = constantValues.RESULTS_FOLDER + constantValues.TEMP_STATS_WORKSPACE + "/" + variableIdentity + str(partitionId) + str(bufferSize) + "/zonalScratch"
    if not os.path.exists(scratchWorkspace): os.makedirs(scratchWorkspace)
    gc.enable()    
    try:
      i=1
      while(i>0):
        tempAddOn = "/zoneTemp" + str(i)
        tempWorkspace = constantValues.RESULTS_FOLDER + constantValues.TEMP_STATS_WORKSPACE + "/" + variableIdentity + str(partitionId) + str(bufferSize) + tempAddOn
        if os.path.exists(tempWorkspace):
          print("scratch workspace " + tempWorkspace + " exists, increasing i value")
          sys.stdout.flush()
          i+=1 
        else:
          print("scratch workspace " + tempWorkspace + "doesn't exist, creating workspace")
          try:
            os.makedirs(tempWorkspace)
            print("sucessfully created " + tempWorkspace)
            i=0
          except:
            i+=1
      tempWorkspace2 = constantValues.RESULTS_FOLDER + constantValues.TEMP_STATS_WORKSPACE + "/" + variableIdentity + str(partitionId) + str(bufferSize) + "/zonalTemp2"
      if not os.path.exists(tempWorkspace): os.makedirs(tempWorkspace) 
      
    except:
      print("couldn't define the tempworkspace for " + variableIdentity)
      sys.stdout.flish()
    try:
      try:
        arcpy.env.scratchworkspace = tempWorkspace2
      except:
        print("couldn't set up the scratchworkspace for " + variableIdentity)
        sys.stdout.flush()
      try:
        arcpy.env.workspace = tempWorkspace
      except:
        print("couldn't define the environmental workspace for " + variableIdentity)
        sys.stdout.flush()
    except:
      print("yo mama")
      sys.stdout.flush()
    try:
      # Check out any necessary licenses
      arcpy.CheckOutExtension("Spatial")
      # Set geoprocessor object property to overwrite existing output, by default
      arcpy.env.overwriteOutput = True
    
      # Local variables...
      inputFeatureZone = inputFile1
      joinField = zoneField
      valueRaster = valueRaster
      joinedFC = outputFile
    
      #Loop through features
      inRows = arcpy.SearchCursor(inputFeatureZone)
      inRow = inRows.next()
    
      #input variables
      strFID = '"FID" ='
      zone = "zone"
    
      cellAssign = "CELL_CENTER"
      priority = "NONE"
      newSelection = "NEW_SELECTION"
      # Make a layer from the feature class
    except:
      print("couldn't define variables")
      sys.stdout.flush()
    try:
      a = arcpy.MakeFeatureLayer_management(inputFeatureZone, zone) 
    except:
        print("couldn't make the feature layer")
        sys.stdout.flush()
            
          
      
    try:
      while inRow:
        try:
          #selct the fid to process
          file = open(constantValues.RESULTS_FOLDER + constantValues.TEMP_STATS_WORKSPACE + "/" + constantValues.TEST_PROGRESS_FILE, 'w')
          file.write("processing variable" + variableIdentity)
          file.close()          
          try:
            a = arcpy.SelectLayerByAttribute_management(zone, newSelection, strFID + str(inRow.FID))
          except:
            print("could not select attribute")
            sys.stdout.flush()
          #create a unique name for the zone
          try:
            uniqueZone = arcpy.CreateUniqueName("zone.shp", arcpy.env.workspace)
            uniqueTable = arcpy.CreateUniqueName("ZSasT", arcpy.env.workspace)
          except:
            print("could not create unique Zone and table")
            sys.stdout.flush()
          #create a temporary feature to use for zonal stats as table
          try:
            arcpy.CopyFeatures_management(zone, uniqueZone)
          except:
            print("could not copy features")
            sys.stdout.flush()
          try:
            outZSaT = ZonalStatisticsAsTable(uniqueZone, joinField, valueRaster, uniqueTable, "DATA", "ALL")
          except Exception as e:
            print("could not run zonal statistics: " + str(e))
            sys.stdout.flush()
          #move to next record.
          try:
            inRow = inRows.next()
          except:
            print("could not retrieve next row")
            sys.stdout.flush()
        except:
          print("could not run the in row while loop")
          sys.stdout.flush()
    except:
      print("failed to run the inRow process, returning to Buffer variables")
      sys.stdout.flush()
    
    #Clear the last selection.
    arcpy.SelectLayerByAttribute_management(zone, "CLEAR_SELECTION")
    #merge the tables so the can be joined to the zone feature class
    mergedTable = arcpy.CreateUniqueName("merge.dbf", arcpy.env.workspace)
    tableList = arcpy.ListTables()
    fullPathTableList = []
    for table in tableList:
        fullPathTableList.append(arcpy.Describe(table).catalogPath)
        
    arcpy.Merge_management(fullPathTableList, mergedTable)
    
    #join them
    arcpy.AddJoin_management(zone, joinField, mergedTable, joinField)
    #save the joined data as a new feature
    arcpy.CopyFeatures_management(zone, joinedFC)
    
    # clean up old output features
    del tableList, zone, outZSaT, uniqueTable, uniqueZone, mergedTable
    print("completed statistics for overlapping zones")
    gc.collect()
    continueVar = 0
  except:
    print("running statistics overlap script failed")
    sys.stdout.flush()
    continueVar = 1
  finally:
    del inRows, inRow
    return continueVar


# end of statistics_for_overlapping_zones.py #