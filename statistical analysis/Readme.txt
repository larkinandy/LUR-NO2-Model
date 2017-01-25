Readme.txt to accompany statistical analysis scripts
# Author: Andrew Larkin
# Created for Perry Hystad, Oregon State University
# Date last modified: January 25, 2017

Statiscal analysis scripts are for creating LUR models in R and evaluating model performance.  

cross_validation.R - used to perform a bootstrap leave 10% out cross-validation.
model_selecdtion.R - perform Lasso variable selecdtion, reduction of selections with 
                     similar buffer sizes, and evaluate selected model performance.
sensitivity_analysis.R - create LUR models based on regional subsets of the continental data.
                         Also create models by fitting residuals of the final global model.