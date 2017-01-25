Readme.txt to accompany Upwind Estimates scripts
# Author: Andrew Larkin
# Created for Perry Hystad, Oregon State University
# Date last modified: January 25, 2017

The python scripts in Upwind Estimates are designed to calcuate the length of major roads upwind from an air monitor on average for a given year
Scripts use a combination of Google Earth Engine and arcpy.  To calculate upwind roads, a buffer file containing road intersections for the buffer area of interest is needed

calculate_upwind_matrix.py - combines road network data with wind direction measurements to estimate length of roads upwind from air monitor
get_wind_vectors_v2.py - download wind direction estimates from GEE for a list of points in lat/long format.

# end of Readme.txt