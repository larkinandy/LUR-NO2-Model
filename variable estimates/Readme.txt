Readme.txt to accompany LUR variable scripts
# Author: Andrew Larkin
# Created for Perry Hystad, Oregon State University
# Date last modified: January 25, 2017

The python scripts in LUR varaible scripts are designed to calcualte point and buffer variable values for a land use regression model.
Multiple values are calculated in parallel for a single buffer distance.  Intermediate products are saved in temporary folders.

buffer_varaiables.py - contains functions used to create buffers and calculate buffer values for both raster and polyline varaible files
constant_values.py - defines contants such as buffer distances, filepaths, and number of simulatenous threads
statistics_for_overlappying_zones.py - creates buffer estimates for a single shapefile and single varaible sequentially, as arcpy does not accurately
calculate values for overlapping buffers within the same shapefile.  This is a modified version of a script provided by NOAA.
run_scripts_v2.py - main python script.  Calls all other scripts to partition the large shapefile in managemenable segments, create buffer and point estimates,
and finally merge partitions back into a single file.

# end of Readme.txt