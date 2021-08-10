## Jericho Lawson
## Summer 2019, 2021
## AFib Functions for Features Program

# Functions that allow the STAGE2features.R code to process the data, 
# create 30 second segments, and add features to the segments.

#### Functions ####

## Utility Functions ##

# Function to determine if value is between two points.
between = function(val, lims){
  return(ifelse(val >= lims[1] & val <= lims[2], T, F))
}

# Function to process each raw atrial fibrillation file for every subject.
file_process = function(number, version){
  if (version == "M"){
    number = paste("0", number, sep = "")
    read.csv(paste("Data/MIT-BIH/pre/", number, ".csv", sep = ""), header = T)
  } else {
    read.csv(paste("Data/2017/pre/", number,  ".csv", sep = ""), header = T)
  }
}


omitter = function(data){
  return(na.omit(data))
}

# Identifies each type that is found in the sample.
type_adder = function(type, counter){
  ind = ifelse(type == "AFIB", 1, ifelse(type == "N", 2, ifelse(type == "AFL", 3, 4)))
  counter[ind] = counter[ind] + 1
  return(counter)
}

state_to_ind = function(type){
  return(ifelse(type == "S", 1, ifelse(type == "Reg", 2, 3)))
}

transition_determiner = function(prev, curr, counter){
  prev_ind = state_to_ind(prev)
  curr_ind = state_to_ind(curr)
  ind = 3 * (prev_ind - 1) + curr_ind
  counter[ind] = counter[ind] + 1
  return(counter)
}


## SAMPLE ENTROPY ##

# finds standardized SampEnt for segment of RR intervals
sampEntStand = function(x, m = 2, r = 0.1 * sdrr, n = length(x)){
  # solve for A and B
  A = sampEntCheck(x, m + 1, r, n)
  B = sampEntCheck(x, m, r, n)
  
  #return(c(A, B)) # debugging purposes
  
  # return standardized sample entropy
  return(-log(A[3] / B[3]))
}

# finds the number of pairs that are lower than threshold for SampEnt
sampEntCheck = function(x, m, r, n){
  # tracks successful differences
  count = 0
  
  # finds max possible differences
  max = choose(n - m + 1, 2)
  
  # pairs each possible m-pack to another m-pack
  for (i in 1:(n + 1 - m - 1)){
    i_pair = x[i:(i + m - 1)]
    for (j in (i + 1):(n + 1 - m)){
      j_pair = x[j:(j + m - 1)]
      count = count + sum(max(abs(i_pair - j_pair)) <= r)
    }
  }
  
  # returns count, maximum, and rate
  return(c(count, max, count / max))
}

## NEC RATE ##

# function to do plus/minus operation
plusminus = function(x, diff){
  if (length(x) != 2){
    NULL
  } else {
    c(x[1] - diff, x[2] + diff)
  }
}

gridplace = function(x, y, int = 0.025){
  xl = plusminus(range(x), (int - (diff(range(x)) %% int)) / 2)
  yl = plusminus(range(y), (int - (diff(range(y)) %% int)) / 2)
  
  # grid lines
  xcuts = seq(xl[1], xl[2], by = int); ycuts = seq(yl[1], yl[2], by = int)
  
  # cell counts
  grid = matrix(0, length(ycuts) - 1, length(xcuts) - 1)
  
  # tally up counts of each cell
  for (p in 1:length(x)){
    coords = c(min(max(1, ceiling((x[p] - xl[1]) / int)), dim(grid)[2]), 
               min(max(1, ceiling((yl[2] - y[p]) / int)), dim(grid)[1]))
    # max added in case point is lower bound for x or upper bound for y

    grid[coords[2], coords[1]] = grid[coords[2], coords[1]] + 1
  }
  return(sum(grid != 0))
}



# Creates several new variables and write data regarding the transition matrix of the AFIB data.

segment_creator = function(subject, list_data, source, chop_left = 0, chop_right = 0, version, int, int_min){
  data = list_data[[subject]]
  # features: transition counter, lengths, differences
  class_count = rep(0, 9)
  
  ind_start = 2 + chop_left
  ind_end = dim(data)[1] - chop_right
  seg_start = ind_start - 1
  
  lengths = c(data$RRLength[seg_start])
  differences = c(data$RRDiff[seg_start])
  
  start = data$Start[seg_start]
  end = data$End[seg_start]
  
  states = c()
  
  if (ind_end - ind_start < 3){ # to run sample entropy properly for 2017 data (n = 5)
    return(NA)
  }
  
  # determines if new data frame is needed
  begin = T
  
  for (i in ind_start:ind_end){
    
    # Adds onto certain parameters and variables.
    end = data$End[i]
    
    # Tracks the amount of each type there is in the sample.
    types_count = type_adder(data$State[i], c(0, 0, 0, 0))
    
    # Tracks the lengths and differences of the RR intervals.
    lengths = c(lengths, data$RRLength[i])
    differences = c(differences, data$RRDiff[i])
    
    # Determines which transition occurs based off the RRclass data.
    class_count = transition_determiner(data$RRClass[i - 1], data$RRClass[i], class_count)
    
    if ((end - start > int & version == "M") | (i == ind_end & version == "C")){
      intervals = i - seg_start + 1
      
      if ((i - seg_start + 1) > int_min) {
        state = ifelse(types_count[1] >= sum(types_count[2:4]), "AFIB", "N") # determines if segment is AFIB or not
        states = c(states, state)
        
        
        heart_rate = intervals * 60 / (end - start);
        props = round(class_count / (sum(class_count)), digits = 8)
        rr_length_var = var(lengths)
        rr_diff_var = var(differences)
        rr_diff_mean = mean(differences)
        
        samp_ent = sampEntStand(x = lengths, m = 2, r = 0.1 * sqrt(rr_length_var))
        necRate = gridplace(x = lengths, y = differences) / intervals
        
        row = c(start, end, (end - start), intervals, heart_rate, 
                props, rr_length_var, rr_diff_var, rr_diff_mean, samp_ent, necRate)
        
        if (begin == T){
          
          rel_freq_info = as.data.frame(rbind("1" = row))
          
          begin = F
        } else {
    
          rel_freq_info = rbind(rel_freq_info, row)
          
        }
      } 
      
      
      start = data$End[i]
      seg_start = i + 1
      types_count = c(0, 0, 0, 0)
      class_count = rep(0, 9)
      lengths = c()
      differences = c()
    }
    
    
    
  }
  # Writes .csv file with info for particular sample.
  if (begin != T & version == "M"){
    rel_freq_info = cbind(states, rel_freq_info)
    class_names = c("State", "Start", "End", "Time", "Intervals", "HeartRate", 
                    "S-S", "S-Reg", "S-L", "Reg-S", "Reg-Reg", "Reg-L", "L-S", 
                    "L-Reg", "L-L", "RRVar", "dRRVar", "dRRMean", "SampEnt", "NECRate")
    colnames(rel_freq_info) = class_names
    
    write.csv(rel_freq_info, paste(source, "seg/", subject, ".csv", sep = ""), row.names = FALSE)
    return(rel_freq_info)
  } else if (begin != T & version == "C") {
    rel_freq_info = cbind(states, rel_freq_info)
    class_names = c("State", "Start", "End", "Time", "Intervals", "HeartRate", 
                    "S-S", "S-Reg", "S-L", "Reg-S", "Reg-Reg", "Reg-L", "L-S", 
                    "L-Reg", "L-L", "RRVar", "dRRVar", "dRRMean", "SampEnt", "NECRate")
    colnames(rel_freq_info) = class_names
    return(rel_freq_info)
  } else {
    return(NA)
  }
  
  
}

combine_to_csv = function(list, source, version, subjects){
  vers = ifelse(version == "M", "MIT-BIH", "2017")
  for (i in 1:length(list)){
    if (length(dim(data2[[i]])) != 0){
      Subject = rep(subjects[i], dim(list[[i]])[1])
      list[[i]] = cbind(Subject, list[[i]])
      if (i == 1){
        to_return = list[[i]]
      } else {
        to_return = rbind(to_return, list[[i]])
      }
    }
  }
  write.csv(to_return, paste(source, "seg/all_seg_data_", vers, ".csv", sep = ""), row.names = FALSE)
}

# Creates several new variables and write data regarding the transition matrix of the AFIB data.
rel_freq_2 = function(num, VERSION){
  # Initializes certain parameters and variables.
  class_count = rep(0, 9); i = 2; start = data[[num]]$Start[1]; end = data[[num]]$End[1]; secs = end - start; intervals = 0
  begin = T
  class_names = c("Start", "End", "Time", "State", "Intervals", "Est. Heart Rate", "S-S", "S-Reg", "S-L", "Reg-S", "Reg-Reg", "Reg-L", "L-S",  
                  "L-Reg", "L-L", "Subject", "R-R Length Variance", "R-R Difference Variance", "R-R Mean Difference")
  
  lengths = c(data[[num]]$RRLength[1]); differences = c(data[[num]]$RRDiff[1])
  
  # Collects dRR interval information.
  drr_all = list()
  
  # Cycles through each heartbeat.
  if (dim(data[[num]])[1] > 1){
    for (i in 2:dim(data[[num]])[1]){
      # Adds onto certain parameters and variables.
      end = data[[num]]$End[i]; secs = secs + data[[num]]$RRLength[i]; intervals = intervals + 1
      
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
      rel_freq_info = rbind(c(start, end, secs, ifelse(data[[num]]$State[1] == "(N", "N", "AFIB"), intervals, intervals * 60 / secs,
                              round(class_count / (sum(class_count)), digits = 8), paste(subjects[num]), 
                              var(lengths), var(differences), mean(differences)))
      drr_all[[num]] = differences
    }else{
      rel_freq_info = rbind(c(start, end, secs, ifelse(data[[num]]$State[1] == "(N", "N", 
                                                     ifelse(data[[num]]$State[1] == "(A", "AFIB",
                                                            ifelse(data[[num]]$State[1] == "(O", "O", "~"))), 
                              intervals, intervals * 60 / secs, round(class_count / (sum(class_count)), digits = 8), paste(subjects[num]), 
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