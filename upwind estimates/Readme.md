# LUR-NO2-Model #

Readme to accompany upwind estimates scripts

**Author:** [Andrew Larkin](https://www.linkedin.com/in/andrew-larkin-525ba3b5/) <br>
**Affiliation:** [Oregon State University, College of Public Health and Human Sciences](https://health.oregonstate.edu/) <br>
**Principal Investigator:** [Perry Hystad](https://health.oregonstate.edu/people/perry-hystad) <br>
**Date Last Modified:** September 23rd, 2018

**Summary** <br>
The python scripts in Upwind Estimates are designed to calcuate the length of major roads upwind from an air monitor on average for a given year. Scripts use a combination of Google Earth Engine and arcpy.  To calculate upwind roads, a buffer file containing road intersections for the buffer area of interest is needed.  

**Processing Scripts** <br>
[**calculate_upwind_matrix.py**](https://github.com/larkinandy/LUR-NO2-Model/blob/master/upwind%20estimates/calculate_upwind_matrix.py) - combines road network data with wind direction measurements to estimate length of roads upwind from air monitor <br>
[**get_wind_vectors_v2.py**](https://github.com/larkinandy/LUR-NO2-Model/blob/master/upwind%20estimates/get_wind_vectors_v2.py) - download wind direction estimates from GEE for a list of points in lat/long format.
