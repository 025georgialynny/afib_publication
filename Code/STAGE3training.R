## Jericho Lawson
## Summer 2019, 2021
## AFib - Resampling and Methods

# Runs machine learning algorithms to classify AFib.

#### SPECIFICATIONS ####################################################################################################

VERSION = "C"   # M for MIT-BIH, C for 2017 CinC dataset
K = 5           # folds of cross-validation (only for 2017 preferably)
SEED = 3        # seed for running machine learning algorithms
DIR = "C:/Users/mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Summer 2021 Testing/" # directory

# Algorithms to run
# LOG = logistic regression
# LDA = linear discriminant analysis
# QDA = quadratic discriminant analysis
# GBM = gradient boosting
# XGB = XGBoost
# LGB = light gradient boosting
# SVM = support-vector machine
# RFO = random forest
# vector of codes below (for TESTS)
TESTS = c("LOG", "LDA", "QDA", "GBM", "XGB", "LGB", "SVM", "RFO")
#TESTS = c("SVM")

# Indices of covariates to use in model
# (Intervals, HeartRate, S.S, S.Reg, S.L, Reg.S, Reg.Reg, Reg.L, 
#      6          7       8     9     10    11      12      13
# L.S, L.Reg, L.L, RRVar, dRRVar, dRRMean, SampEnt, CoSEn, NECRate)
#  14   15    16    17      18      19       20       21     22
COVARIATES = c(8:16, 17, 18, 22)

# Index of response vector (State)
RESPONSE = 2

#### DATA GATHERING & PREPARATION ######################################################################################

# Sets directory; gathers file with functions to run program; opens libraries needs to run models.
setwd(DIR)
source("Code/STAGE3functions.R")
sapply(TESTS, library_selector)

# Retrieves all segment information for specific dataset.
if (VERSION == "M"){ # "Data/MIT-BIH/seg" for MIT-BIH, "Data/2017/seg" for 2017 CinC dataset
  subpath = "Data/MIT-BIH/seg/"
  samp_data = read.csv(paste(subpath, "all_seg_data_MIT-BIH.csv", sep = ""))
} else {
  subpath = "Data/2017/seg/"
  samp_data = read.csv(paste(subpath, "all_seg_data_2017.csv", sep = ""))
}

# Converts state column from strings to factors.
samp_data$State = factor(samp_data$State, levels = c("N", "AFIB"))

# Creates formula for models.
form = formula_maker(COVARIATES, colnames(samp_data))

#### CROSS-VALIDATION AND TESTING ######################################################################################

# Randomly assigns number to each observation for k-fold cross validation.
groups = grouper(samp_data, VERSION, SEED, K)

# For each algorithm, we use the testing function to run k-fold cross validation on the specified model
res = lapply(TESTS, testing, data = samp_data, rand = groups, formula = form, index_c = COVARIATES, index_r = RESPONSE, 
           thresh = 0.5, seed = SEED)

# Add names to the list of results.
names(res) = TESTS

# Calculate metrics indicating the success of the models and algorithms.
metrics = lapply(res, metrics_f)
metrics

########################################################################################################################
