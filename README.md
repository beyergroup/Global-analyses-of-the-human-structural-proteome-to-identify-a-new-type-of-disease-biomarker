# Global analyses of the human structural proteome to identify a new type of disease biomarker

Scripts to reproduce the main analysis of "Global analyses of the human structural proteome to identify a new type of disease biomarker"



## Main analysis

**StructuralVariabiliyCSF.R**: Variability of the healthy human CSF structural proteome

**LinearModels.R**: Structural proteomic changes in CSF of PD patients

**RegularizedRegressionModels.R**: CSF structural peptides classify healthy vs PD groups

## Requirements

### Data requirements
Please download supplementary file 12, all data you need for running the analysis is provided there.


### System requirements

The code runs on Linux operating systems (Fedora 25) and has not been tested on other versions.

The code was developed on R-version 4.0.3 including following packages:
```
devtools 2.4.3
glmnet 4.1-1
lawstat 3.4
multimode 1.5
openxlsx 4.2.5
plyr 1.8.7 
```
