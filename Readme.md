![GitHub Logo](/GlobalLUR.jpg)

**Author:** [Andrew Larkin](https://www.linkedin.com/in/andrew-larkin-525ba3b5/) <br>
**Affiliation:** [Oregon State University, College of Public Health and Human Sciences](https://health.oregonstate.edu/) <br>
**Principal Investigator:** [Perry Hystad](https://health.oregonstate.edu/people/perry-hystad) <br>
**Date last modified:** September 23rd, 2018

**Summary:** <br>
This github repository contains the scripts used to create a global NO2 land use regression model. Land use estimates were derived using python and ArcGIS.  Variable selection and model development were developed using R Studio.

**Repository Structure:** <br>
Files are divided into three folders, with each folder corresponding to a unique stage of model development.

- **[variable estimates](https://github.com/larkinandy/LUR-NO2-Model/tree/master/variable%20estimates)** - scripts in this folder were used to derive NO2 estiamtes at air monitor locations using parallel processing. <br>
- **[upwind estimates](https://github.com/larkinandy/LUR-NO2-Model/tree/master/upwind%20estimates)** - scripts in this folder were used to estimate length of road upwind from air monitor locations.
- **[statistical analysis](https://github.com/larkinandy/LUR-NO2-Model/tree/master/statistical%20analysis)** - scripts in this folder were used to perform lasso variable selection, model evaluation, and sensitivity analysis. <br>

**External Links**
- **Publication** - https://pubs.acs.org/doi/abs/10.1021/acs.est.7b01148
- **Model estimates and underlying datasources** - https://health.oregonstate.edu/labs/spatial-health/resources-equipment
