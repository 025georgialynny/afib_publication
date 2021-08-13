## Jericho Lawson
## Summer 2019, 2021
## AFib - Resampling and Methods

# Runs machine learning algorithms to classify AFib.

#### SPECIFICATIONS ####################################################################################################

VERSION = "M" # M for MIT-BIH, C for 2017 CinC dataset
K = 5 # folds of cross-validation (only for 2017 preferably)
SEED = 3 # seed for running machine learning algorithms
DIR = "C:/Users/mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Summer 2021 Testing/"

# LOG = logistic regression, LDA = linear discriminant analysis
# QDA = quadratic discriminant analysis, GBM = gradient boosting
# XGB = XGBoost, LGB = light gradient boosting
TESTS = c("LOG", "LDA", "QDA")#, "GBM") 

# (Intervals, HeartRate, S.S, S.Reg, S.L, Reg.S, Reg.Reg, Reg.L, 
#      6          7       8     9     10    11      12      13
# L.S, L.Reg, L.L, RRVar, dRRVar, dRRMean, SampEnt, NECRate)
#  14   15    16    17      18      19       20       21
COVARIATES = 21

#### DATA COLLECTION ###################################################################################################

# Sets directory; gathers file with functions to run program.
setwd(DIR)
source("Code/STAGE3functions.R")

# Retrieves all R-R interval information for all subjects.
if (VERSION == "M"){ # "Data/MIT-BIH/seg" for MIT-BIH, "Data/2017/seg" for 2017 CinC dataset
  subpath = "Data/MIT-BIH/seg/"
  samp_data = read.csv(paste(subpath, "all_seg_data_MIT-BIH.csv", sep = ""))
} else {
  subpath = "Data/2017/seg/"
  samp_data = read.csv(paste(subpath, "all_seg_data_2017.csv", sep = ""))
}

#### CROSS-VALIDATION AND TESTING ######################################################################################

samp_data$State = factor(samp_data$State, levels = c("N", "AFIB"))
groups = grouper(samp_data, VERSION, SEED, K)
form = formula_maker(COVARIATES, colnames(samp_data))

res = lapply(TESTS, testing, data = samp_data, rand = groups, formula = form, 
           thresh = 0.5, seed = SEED)

metrics(res[[2]])

########################################################################################################################
