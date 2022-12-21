# data_200_final_project
Tufts University

# Contents:
## data
This folder contains all the files needed to replicate the analysis for the series of Forward Stepwise Regression models. The `diet.xpt` file is the result of the data processed version of the [DR1IFF_I](https://wwwn.cdc.gov/nchs/nhanes/search/datapage.aspx?Component=Dietary&CycleBeginYear=2015) file available for download on the NHANES website. That file was too large to be uploaded to GitHub.

### extra_data
This subdirectory contains additional information to the potential predictor variables for hypertension. This contains the unprocessed medication data as well as the file containing our response variables (systolic and diastolic blood pressure measurements of each individual).

## forward_selection_final.R
I annotated the R Script used to run the Forward Stepwise Regression models as well as performed any data cleaning needed prior to running the model.

## data_visualization_plots_R.R
I produced Figures 1 and 2 of our report using the annotated code that I provided in this file.

## saved_plots
This directory contains the saved plots (Figures 1 and 2) encoded in the `data_visualization_plots_R.R` file.
