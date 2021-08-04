## Jericho Lawson
## Summer 2019
## Atrial Fib
## Used to create features for AFib classification.

#### VERSION & GROUP #### (for user)

# Pick specific version and group; group parameters can be changed in the GLOBALS section. Groups allow you to change the
# thresholds for cleaning some out of the outlier data.
VERSION = "MIT"
#VERSION = "single_2017"
#VERSION = "multi_2017"
GROUP = ""
#GROUP = "_1"
#GROUP = "_2"

# Decide whether you want to add covariates from the Kolmogorov-Smirnov tests (for MIT-BIH data only).
#KS = "yes"
KS = "no"

###################################################################################################

#### Libraries and Globals ####

library(ggplot2); library(gridExtra); library(wesanderson); library(reshape2); library(ggpubr)
library(tidyr)

#### GLOBALS ####
if (GROUP == ""){
  W1 = 0.75; W2 = 0.25; INTERVAL = 30; TYPE_THRESH = "n"; RUN_THRESH = c(0.85, 1.15) # DO NOT CHANGE
}else if (GROUP == "_1"){
  W1 = 0.67; W2 = 0.33; INTERVAL = 30; TYPE_THRESH = "n"; RUN_THRESH = c(0.9, 1.1)
}else{
  W1 = 0.75; W2 = 0.25; INTERVAL = 30; TYPE_THRESH = "n"; RUN_THRESH = c(0.85, 1.15)
}

#### Data Gathering ####

source("C:/Users/Mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Summer 2021 Testing/Code/STAGE2atrial_fib_functs.R")

# Retrieves all R-R interval information for all subjects.
if (VERSION == "MIT"){
  subjects = c('04015', '04043', '04048', '04126', '04746', '04908', '04936', '05091', 
               '05121', '05261', '06426', '06453', '06995', '07162', '07859', '07879', 
               '07910', '08215', '08219', '08378', '08405', '08434', '08455')
}else{
  subjects = read.csv("C:/Users/mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Testing/training2017/REFERENCE_2.csv", header = F)$V2[-1]
}
data = lapply(subjects, file_process)

#### Data Processing ####

# Determines running mean and class for each observation for each sample in data.
for (item in 1:length(data)){
  RR_mean = rep(NA, dim(data[[item]])[1]); RR_class = rep(NA, dim(data[[item]])[1]); diffs = 0; 
  RR_diff = rep(NA, dim(data[[item]])[1])
  for (obs in 1:dim(data[[item]])[1]){
    
    # Determines if observation is an outlier.
    res = outlier(TYPE_THRESH, 2, data[[item]]$RRLength[obs], data[[item]]$RRLength)
    
    # Determines a class and running mean for the observation if it is not an outlier.
    if (res == T){
      # Used for first observation in data.
      if (obs - diffs == 1){
        RR_mean[obs] = data[[item]]$RRLength[obs]
        RR_diff[obs] = 0
        diffs = 0
      }else{
        # Used if there are any holes in the data due to the outliers.
        RR_mean[obs] = W1 * RR_mean[obs - diffs - 1] + W2 * data[[item]]$RRLength[obs]
        RR_diff[obs] = data[[item]]$RRLength[obs] - data[[item]]$RRLength[obs - diffs - 1]
        if (diffs != 0){
          diffs = 0
        }
      }
      # Determine class based off how far away the length is from the running mean.
      RR_class[obs] = ifelse(data[[item]]$RRLength[obs] <= (RUN_THRESH[1] * RR_mean[obs]), "S",
                             ifelse(data[[item]]$RRLength[obs] <= (RUN_THRESH[2] * RR_mean[obs]), 
                                    "Reg", "L"))
    }else{
      diffs = diffs + 1
    }
  }
  # Assigns a running mean and class to the observation. (alter as needed)
  data[[item]]$RRMean = RR_mean; data[[item]]$RRClass = RR_class; data[[item]]$RRDiff = RR_diff
  if (VERSION == "MIT"){
    write.csv(data[[item]], paste("C:/Users/mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Summer 2021 Testing/Data/", subjects[item], VERSION, 
                                  "expanded", GROUP, ".csv", sep = ""))
  }else if (VERSION == "single_2017"){
    write.csv(data[[item]], paste("C:/Users/mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Summer 2021 Testing/Data/", subjects[item], VERSION,
                                  "expanded", GROUP, ".csv", sep = ""))
  }else{
    write.csv(data[[item]], paste("C:/Users/mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Summer 2021 Testing/Data/", subjects[item], VERSION,
                                  "multi_expanded", GROUP, ".csv", sep = ""))
  }
}

# Omits any NA data.
for (item in 1:length(data)){
  data[[item]] = na.omit(data[[item]])
}

# Compile the data together and write a new .csv file.
if (VERSION == "MIT"){
  data_2_pre = lapply(1:length(subjects), rel_freq, VERSION = VERSION)
  
  # Used for creating covariates from the Kolmogorov-Smirnov Tests.
  ks_list = list()
  
  # Used to translate characters into numerics.
  for (item in 1:length(data_2_pre)){
    
    ks_list[[item]] = list(); count = 1
    
    for (vector in 1:length(data_2_pre[[item]][[2]])){
      ks_list[[item]][[count]] = data_2_pre[[item]][[2]][[vector]]; count = count + 1
    }
    
    data_2 = as.data.frame(data_2_pre[[item]][[1]])
    data_2$Subject = subjects[item]
    if (item == 1){
      all_data_2 = rbind(data_2)
    }else{
      all_data_2 = rbind(all_data_2, data_2)
    }
  }
  write.csv(all_data_2, paste("C:/Users/mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Summer 2021 Testing/Data/all_trans_data", VERSION, GROUP, ".csv", sep = ""))
  
  # Creation of Kolmogorov-Smirnov test matrices for all users.
  if (KS == "yes"){
    ks_matrices = list()
    for (m in 1:length(data_2_pre)){
      ks_mat = matrix(NA, length(data_2_pre[[m]][[2]]), length(data_2_pre[[m]][[2]]))
      j = 1
      k = 1
      for (i in 1:length(ks_list[[m]])){
        while (k <= length(ks_list[[m]])){
          if (i != k){
            ks_mat[k, i] = ks.test(ks_list[[m]][[i]], ks_list[[m]][[k]])$statistic
          }
          ks_mat[i, k] = ks.test(ks_list[[m]][[i]], ks_list[[m]][[k]])$statistic
          k = k + 1
        }
        j = j + 1
        k = j
      }
      write.csv(ks_mat, paste("C:/Users/mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Summer 2021 Testing/Data/", subjects[m], "ks.csv", sep = ""))
    }
  }
  
}else{
  data_2_pre = lapply(1:length(data), rel_freq_2, VERSION = VERSION)
  
  # Used for creating covariates from the Kolmogorov-Smirnov Tests. (Can be coded to add the actual K-S tests in this dataset.)
  ks_list = list(); count = 1
  
  for (item in 1:length(data_2_pre)){
    if(length(data_2_pre[[1]][[item]]) != 1){
      data_2 = as.data.frame(data_2_pre[[item]][[1]])
      
      ks_list[[count]] = data_2_pre[[item]][[2]]; count = count + 1
      
      if (item == 1){
        all_data_2 = rbind(data_2)
      }else{
        all_data_2 = rbind(all_data_2, data_2)
      }
    }
  }
  if (VERSION == "single_2017"){
    write.csv(all_data_2, paste("C:/Users/mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Summer 2021 Testing/Data/all_tran_mat_2017", VERSION, GROUP, ".csv", sep = ""))
  }else{
    write.csv(all_data_2, paste("C:/Users/mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Summer 2021 Testing/Data/all_tran_mat_2017_multi", VERSION, GROUP, ".csv", sep = ""))
  }
}

###################################################################################################