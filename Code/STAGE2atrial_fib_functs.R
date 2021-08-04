## Jericho Lawson
## Summer 2019
## Atrial Fib Functions for Atrial Fib Program

#### Functions ####

VERSION = "MIT"
#VERSION = "single_2017"
#VERSION = "multi_2017"

# Function to process each raw atrial fibrillation file for every subject.
file_process = function(number){
  if (VERSION == "MIT"){
    read.csv(paste("C:/Users/mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Summer 2021 Testing/Data/", number, 
                 ".csv", sep = ""), header = T)
  }else if (VERSION == "single_2017"){
    read.csv(paste("C:/Users/mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Data/Training_Set/", number, 
                   ".csv", sep = ""), header = T)
  }else{
    read.csv(paste("C:/Users/mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Data/Training_Set/", number, 
                 "_multi.csv", sep = ""), header = T)
  }
}

# Determines if a heartbeat is an outlier or not based off a certain threshold.
outlier = function(np, thresh, length, length_data){
  return(ifelse(np == "n", ifelse(length >= thresh, F, T), ifelse(length >= quantile(length_data, thresh), F, T)))
}

# Identifies each type that is found in the sample.
types_funct = function(type, types_count){
  if (type == "(AFIB"){
    types_count[1] = types_count[1] + 1
  }else if (type == "(N"){
    types_count[2] = types_count[2] + 1
  }else if (type == "(AFL"){
    types_count[3] = types_count[3] + 1
  }else{
    types_count[4] = types_count[4] + 1
  }
  return(types_count)
}

# Creates several new variables and write data regarding the transition matrix of the AFIB data.
rel_freq = function(num, VERSION){
  # Initializes certain parameters and variables.
  class_count = rep(0, 9); i = 2; start = data[[num]]$Start[1]; end = data[[num]]$End[1]; secs = end - start; heartbeats = 0
  count = 1
  begin = T; types_count = types_funct(data[[num]]$State[1], c(0, 0, 0, 0))
  class_names = c("Start", "End", "Time", "State", "Heartbeats", "Est. Heart Rate", "S-S", "S-Reg", "S-L", "Reg-S", "Reg-Reg", "Reg-L", "L-S", 
                  "L-Reg", "L-L", "R-R Length Variance", "R-R Difference Variance", "R-R Mean Difference")
  
  lengths = c(data[[num]]$RRLength[1]); differences = c(data[[num]]$RRDiff[1])
  
  # Collects dRR interval information.
  drr_all = list()
  
  # Cycles through each heartbeat.
  while (i < dim(data[[num]])[1]){
    # Adds onto certain parameters and variables.
    end = data[[num]]$End[i]; secs = secs + data[[num]]$RRLength[i]; heartbeats = heartbeats + 1
    
    # Tracks the amount of each type there is in the sample.
    types_count = types_funct(data[[num]]$State[i], types_count)
    
    # Tracks the lengths and differences of the RR intervals.
    lengths = c(lengths, data[[num]]$RRLength[i]); differences = c(differences, data[[num]]$RRDiff[i])
    
    # Determines which transition occurs based off the RRclass data.
    if (data[[num]]$RRClass[i - 1] == "S"){
      if (data[[num]]$RRClass[i] == "S"){
        class_count[1] = class_count[1] + 1
      }else if (data[[num]]$RRClass[i] == "Reg"){
        class_count[2] = class_count[2] + 1
      }else{
        class_count[3] = class_count[3] + 1
      }
    }else if (data[[num]]$RRClass[i - 1] == "Reg"){
      if (data[[num]]$RRClass[i] == "S"){
        class_count[4] = class_count[4] + 1
      }else if (data[[num]]$RRClass[i] == "Reg"){
        class_count[5] = class_count[5] + 1
      }else{
        class_count[6] = class_count[6] + 1
      }
    }else{
      if (data[[num]]$RRClass[i] == "S"){
        class_count[7] = class_count[7] + 1
      }else if (data[[num]]$RRClass[i] == "Reg"){
        class_count[8] = class_count[8] + 1
      }else{
        class_count[9] = class_count[9] + 1
      }
    }
    
    # Places data point into new data file.
    if (secs > INTERVAL){
      state = ifelse(types_count[1] >= sum(types_count[2:4]), "AFIB", "N")
      #state = ifelse(types_count[1] == sum(types_count), "AFIB", ifelse(sum(types_count[2:4]) == sum(types_count), "N", "F"))
      if (state != "F"){
        if (begin == T){
          rel_freq_info = rbind(c(start, data[[num]]$End[i], secs, state, heartbeats, heartbeats * 60 / secs,
                                  round(class_count / (sum(class_count)), digits = 8), var(lengths), var(differences), mean(differences)))
          drr_all[[count]] = differences
          count = count + 1
          begin = F
        }else{
          rel_freq_info = rbind(rel_freq_info, c(start, data[[num]]$End[i], secs, state, heartbeats, heartbeats * 60 / secs,
                                                 round(class_count / (sum(class_count)), digits = 8), 
                                                 var(lengths), var(differences), mean(differences)))
          drr_all[[count]] = differences
          count = count + 1
        }
      }
      start = data[[num]]$End[i]; secs = 0; types_count = c(0, 0, 0, 0); class_count = rep(0, 9); heartbeats = 0
      lengths = c(); differences = c()
    }
    
    # Iterates the index of the sample.
    i = i + 1
  }
  
  # Writes .csv file with info for particular sample.
  if (begin != T){
    colnames(rel_freq_info) = class_names
    write.csv(rel_freq_info, paste("C:/Users/mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Data/", subjects[num], VERSION,
                                 "tran_mat.csv", sep = ""))
    return(list(rel_freq_info, drr_all))
  }else{
    return(NA)
  }
}

# Creates several new variables and write data regarding the transition matrix of the AFIB data.
rel_freq_2 = function(num, VERSION){
  # Initializes certain parameters and variables.
  class_count = rep(0, 9); i = 2; start = data[[num]]$Start[1]; end = data[[num]]$End[1]; secs = end - start; heartbeats = 0
  begin = T
  class_names = c("Start", "End", "Time", "State", "Heartbeats", "Est. Heart Rate", "S-S", "S-Reg", "S-L", "Reg-S", "Reg-Reg", "Reg-L", "L-S",  
                  "L-Reg", "L-L", "Subject", "R-R Length Variance", "R-R Difference Variance", "R-R Mean Difference")
  
  lengths = c(data[[num]]$RRLength[1]); differences = c(data[[num]]$RRDiff[1])
  
  # Collects dRR interval information.
  drr_all = list()
  
  # Cycles through each heartbeat.
  if (dim(data[[num]])[1] > 1){
    for (i in 2:dim(data[[num]])[1]){
      # Adds onto certain parameters and variables.
      end = data[[num]]$End[i]; secs = secs + data[[num]]$RRLength[i]; heartbeats = heartbeats + 1
      
      # Tracks the lengths and differences of the RR intervals.
      lengths = c(lengths, data[[num]]$RRLength[i]); differences = c(differences, data[[num]]$RRDiff[i])
      
      # Determines which transition occurs based off the RRclass data.
      if (data[[num]]$RRClass[i - 1] == "S"){
        if (data[[num]]$RRClass[i] == "S"){
          class_count[1] = class_count[1] + 1
        }else if (data[[num]]$RRClass[i] == "Reg"){
          class_count[2] = class_count[2] + 1
        }else{
          class_count[3] = class_count[3] + 1
        }
      }else if (data[[num]]$RRClass[i - 1] == "Reg"){
        if (data[[num]]$RRClass[i] == "S"){
          class_count[4] = class_count[4] + 1
        }else if (data[[num]]$RRClass[i] == "Reg"){
          class_count[5] = class_count[5] + 1
        }else{
          class_count[6] = class_count[6] + 1
        }
      }else{
        if (data[[num]]$RRClass[i] == "S"){
          class_count[7] = class_count[7] + 1
        }else if (data[[num]]$RRClass[i] == "Reg"){
          class_count[8] = class_count[8] + 1
        }else{
          class_count[9] = class_count[9] + 1
        }
      }
      begin = F
    }
    # Make changes below if necessary. (1st for binary, 2nd for multiple classes)
    if (VERSION == "single_2017"){
      rel_freq_info = rbind(c(start, end, secs, ifelse(data[[num]]$State[1] == "(N", "N", "AFIB"), heartbeats, heartbeats * 60 / secs,
                              round(class_count / (sum(class_count)), digits = 8), paste(subjects[num]), 
                              var(lengths), var(differences), mean(differences)))
      drr_all[[num]] = differences
    }else{
      rel_freq_info = rbind(c(start, end, secs, ifelse(data[[num]]$State[1] == "(N", "N", 
                                                     ifelse(data[[num]]$State[1] == "(A", "AFIB",
                                                            ifelse(data[[num]]$State[1] == "(O", "O", "~"))), 
                              heartbeats, heartbeats * 60 / secs, round(class_count / (sum(class_count)), digits = 8), paste(subjects[num]), 
                              var(lengths), var(differences), mean(differences)))
      drr_all[[num]] = differences
    }
  }
  # Writes .csv file with info for particular sample.
  if (begin == F){
    colnames(rel_freq_info) = class_names
    return(list(rel_freq_info, drr_all))
  }else{
    return(1)
  }
}

# Gives proportion of sample that is AFIB and non-AFIB.
prop_funct = function(num){
  vect = c(subjects[num], sum(ifelse(data[[num]]$State == "(AFIB", "AFIB", "N") == "AFIB") / length(data[[num]]$State), 
           sum(ifelse(data[[num]]$State == "(AFIB", "AFIB", "N") != "AFIB") / length(data[[num]]$State), length(data[[num]]$State), 
           summary(data[[num]]$RRLength))
  return(vect)
}

###################################################################################################