######### constant_values.py ########
# header file to store all constants used in processing buffer variables for LUR models
# author: Andrew Larkin
# created for Perry Hsytad, Oregon State University
# date last modified: 01/25/2017


############## define settings and variables #########
KEYWORD = "Partition"
ZONE_KEYWORD = 'z' # character added to filepath to indicate which zone air monitor partition is located in
PARTITION_KEYWORD = "i" # character added to filepath to indicate the partition number of an air monitor partition
RASTER_TYPE = 0 # unique identifier for files defiend as rasters
POLYLINE_TYPE = 1 # unique identifier for files defined as polylines
POINT_TYPE = 3 # unique identifier for files defined as points that don't require buffer analysis
PARALLEL_PROCESSING = 4 # number of simulatenous threads to spawn
POINT_BUFFER_TYPE = 4 # unique identifier for files defined as points that do require buffer analysis
BUFFER_EXTENSION = "/buffers/" # subfolder containing buffer shapefiles
AIRMONITOR_ID = "FID"
TABLE_ID = "ORIG_FID"
BUFFER_ID = "FID_buffer"
LENGTH_COMMAND = "!shape.length@kilometers!" # command for calculating road length
CARBON_COMMAND = "!SUM_carbon! * 1" # command for calculating power plant CO2 buffer estiamtes
TEMP_STATS_WORKSPACE = "tempStats" # temporary workspace - with multiprocessing, each thread needs a unique workspace
RATER_PROCESS_WAIT_TIME = 60*5 # amount of time to wait between checks to see if parallel threads are still running
TEST_PROGRESS_FILE = "test_progress.txt" # test print statement used to check if threads are running or have stalled
MONITOR_FILE= "nycFishp2.shp" # input air monitor shapefile to process
INPUT_FOLDER = "insert input folder here" # input folder containing both the air monitor file and variable files
RESULTS_FOLDER = "insert results folder here" # folder to store results in, including temporary air monitor partitions
ZONE_DEFINITIONS = "zoneDef.shp" # shapefile used to categorize air monitor locations into specific zones
MOSAIC_RASTER_LIST = [] 
MOSAIC_RASTER_LIST = ["trCov"] # variable files in raster mosaic format that require buffer analysis
POLLYLINE_MOSAIC_LIST = ["mjRds","miRds"] # variable files in mosaic polyline format
POINT_BUFFER_LIST = ["plants.shp"] # variable files in point format that need buffer analysis
POINT_MOSAIC_LIST = [] # variable files in raster mosaic format that do not require buffer analysis
POINT_LIST = ["satNO2.tif"] # variable files in raster format that do not require buffer analysis
RASTER_LIST = ["fireAvg1.tif","popDens4.tif","imp2010.tif","water_percent.tif","NDVI_2010v21.tif"] # variable files in raster format that do require buffer analysis
POLYLINE_LIST = [] # variable files in polyline format that require buffer analysis
COAST_BOUNDARY = "coastline.shp" # file defining coastline
BUFFER_DISTANCE = [100, 200, 300, 400, 500, 600, 700, 800, 1000, 1200, 1500, 2000, 2500, 3000, 3500, 4000, 5000, 6000, 7000, 8000, 10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000, 50000]
PARTITION_SIZE = 20 # size for eaach air monitor partition.  Larger partitions and buffer sizes require computers with more RAM.

####### end of define settings and variables #########



### end of constant_values.py ###