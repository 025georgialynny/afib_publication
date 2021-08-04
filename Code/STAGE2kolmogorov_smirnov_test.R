## Jericho Lawson
## July 12th, 2019
## Kolmogorov-Smirnov Tests

#### Libraries and Globals ####

library(ggplot2); library(gridExtra); library(reshape2); library(ggpubr)
library(tidyr)

#### VERSION ####
VERSION = "MIT"
#VERSION = "single_2017"
#VERSION = "multi_2017"
#GROUP = ""
#GROUP = "_1"
GROUP = "_2"

#### GLOBALS ####
if (GROUP == ""){
  W1 = 0.75; W2 = 0.25; INTERVAL = 30; TYPE_THRESH = "n"; RUN_THRESH = c(0.85, 1.15)
}else if (GROUP == "_1"){
  W1 = 0.67; W2 = 0.33; INTERVAL = 30; TYPE_THRESH = "n"; RUN_THRESH = c(0.9, 1.1)
}else{
  W1 = 0.75; W2 = 0.25; INTERVAL = 30; TYPE_THRESH = "n"; RUN_THRESH = c(0.85, 1.15)
}

#### Data Gathering ####

source("/rhome/jlaws011/AFib/STAGE2atrial_fib_functs_hpc.R")

if (VERSION == "MIT"){
  subjects = c('04015', '04043', '04048', '04126', '04746', '04908', '04936', '05091', 
               '05121', '05261', '06426', '06453', '06995', '07162', '07859', '07879', 
               '07910', '08215', '08219', '08378', '08405', '08434', '08455')
  data = lapply(subjects, file_process)
}else{
  subjects = read.csv("/rhome/jlaws011/AFib/REFERENCE_2.csv", header = F)$V2[-1]
  if (VERSION == "single_2017"){
    data = lapply(subjects, file_process_2)
  }else{
    data = lapply(subjects, file_process_3)
  }
}

#### Data Processing ####

# Determines running mean and class for each observation for each sample in data.
for (item in 1:length(data)){
  RR_mean = rep(NA, dim(data[[item]])[1]); RR_class = rep(NA, dim(data[[item]])[1]); diffs = 0; 
  RR_diff = rep(NA, dim(data[[item]])[1])
  for (obs in 1:dim(data[[item]])[1]){
    
    # Determines if observation is an outlier.
    res = outlier(TYPE_THRESH, 1.5, data[[item]]$RRLength[obs], data[[item]]$RRLength)
    
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
}

# Omits any NA data.
for (item in 1:length(data)){
  data[[item]] = na.omit(data[[item]])
}

# Compile the data together and write a new .csv file.
if (VERSION == "MIT"){
  data_2_pre = lapply(1:length(subjects), rel_freq, VERSION = VERSION)
  
  ks_list = list(); count = 1
  
  # Used to translate characters into numerics.
  for (item in 1:length(data_2_pre)){
    
    for (vector in 1:length(data_2_pre[[item]][[2]])){
      ks_list[[count]] = data_2_pre[[item]][[2]][[vector]]; count = count + 1
    }
    
    data_2 = as.data.frame(data_2_pre[[item]][[1]])
    data_2$Subject = subjects[item]
    if (item == 1){
      all_data_2 = rbind(data_2)
    }else{
      all_data_2 = rbind(all_data_2, data_2)
    }
  }
  #write.csv(all_data_2, paste("C:/Users/mario/Documents/UNCW 2019/Atrial_Fibrillation/Data/all_trans_data", VERSION, GROUP, ".csv", sep = ""))
  
  
  ##New
  ks_matrices = list()
  for (m in 1:length(data_2_pre)){
    ks_mat = matrix(NA, length(data_2_pre[[m]]), length(data_2_pre[[m]]))
    j = 1
    k = 1
    for (i in 1:length(ks_list)){
      while (k <= length(ks_list)){
        if (i != k){
          ks_mat[k, i] = ks.test(ks_list[[i]], ks_list[[k]])$statistic
        }
        ks_mat[i, k] = ks.test(ks_list[[i]], ks_list[[k]])$statistic
        k = k + 1
      }
      j = j + 1
      k = j
    }
    write.csv(ks_mat, paste("/rhome/jlaws011/AFib/KS_DATA/", subjects[m], "ks.csv", sep = ""))
  }
  
    
  
}else{
  data_2_pre = lapply(1:length(data), rel_freq_2, VERSION = VERSION)
  
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
    #write.csv(all_data_2, paste("C:/Users/mario/Documents/UNCW 2019/Atrial_Fibrillation/Data/Training_Set/all_tran_mat_2017", VERSION, GROUP, ".csv", sep = ""))
  }else{
    #write.csv(all_data_2, paste("C:/Users/mario/Documents/UNCW 2019/Atrial_Fibrillation/Data/Training_Set/all_tran_mat_2017_multi", VERSION, GROUP, ".csv", sep = ""))
  }
}

write.csv(ks_mat, "/rhome/jlaws011/AFib/ks_mat.csv", sep = "")



