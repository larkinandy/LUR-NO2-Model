# LUR-NO2-Model #

Readme.txt to accompany LUR variable scripts

**Author:** [Andrew Larkin](https://www.linkedin.com/in/andrew-larkin-525ba3b5/) <br>
**Affiliation:** [Oregon State University, College of Public Health and Human Sciences](https://health.oregonstate.edu/) <br>
**Principal Investigator:** [Perry Hystad](https://health.oregonstate.edu/people/perry-hystad) <br>
**Date last modified:** September 23rd, 2018

**Summary** <br>
The python scripts in LUR varaible scripts are designed to calcualte point and buffer variable values for a land use regression model. Multiple values are calculated in parallel for a single buffer distance.  Intermediate products are saved in temporary folders. 

**Processing Scripts** <br>
- [**buffer_varaiables.py**](https://github.com/larkinandy/LUR-NO2-Model/blob/master/variable%20estimates/buffer_variables.py) - contains functions used to create buffers and calculate buffer values for both raster and polyline varaible files
- [**constant_values.py**](https://github.com/larkinandy/LUR-NO2-Model/blob/master/variable%20estimates/constant_values.py) - defines contants such as buffer distances, filepaths, and number of simulatenous threads
- [**statistics_for_overlappying_zones.py**](https://github.com/larkinandy/LUR-NO2-Model/blob/master/variable%20estimates/statistics_for_overlapping_zones.py) - creates buffer estimates for a single shapefile and single varaible sequentially, as arcpy (as tested in version 10.3.1) does not accurately calculate values for overlapping buffers within the same shapefile. This is a modified version of a script provided by NOAA.
- [**run_scripts_v2.py**](https://github.com/larkinandy/LUR-NO2-Model/blob/master/variable%20estimates/run_scripts_v2.py) - main python script.  Calls all other scripts to partition the large shapefile in managemenable segments, create buffer and point estimates, and finally merge partitions back into a single file.
