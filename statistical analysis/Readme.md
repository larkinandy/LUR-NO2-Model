# LUR-NO2-Model #

Readme to accompany statistical analysis scripts

**Author:** [Andrew Larkin](https://www.linkedin.com/in/andrew-larkin-525ba3b5/) <br>
**Affiliation:**[Oregon State University, College of Public Health and Human Sciences](https://health.oregonstate.edu/) <br>
**Principal Investigator:**[Perry Hystad](https://health.oregonstate.edu/people/perry-hystad) <br>
**Date Last Modified:** September 23rd, 2018

**Summary:** <br>
Scripts in this folder are for creating and evaluating the near surface NO2 LUR model based on variables already derived from scripts in other folders in this github.  

**Processing Scripts** <br>
- [**cross_validation.R**](https://github.com/larkinandy/LUR-NO2-Model/blob/master/statistical%20analysis/cross_validation.R) - used to perform a bootstrap leave 10% out cross-validation.
- [**model_selection.R**](https://github.com/larkinandy/LUR-NO2-Model/blob/master/statistical%20analysis/model_selection.R) - perform Lasso variable selecdtion, reduction of selections with similar buffer sizes, and evaluate selected model performance.
- [**sensitivity_analysis.R**](https://github.com/larkinandy/LUR-NO2-Model/blob/master/statistical%20analysis/sensitivity_analysis.R) - create LUR models based on regional subsets of the continental data. Also create models by fitting residuals of the final global model.
