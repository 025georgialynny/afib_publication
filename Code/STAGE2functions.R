## Jericho Lawson
## Summer 2019, 2021
## AFib Functions for Features Program

# Functions that allow the STAGE2features.R code to process the data, 
# create 30 second segments, and add features to the segments.

#### Helper Functions ##########################################################

# Function to add and subtract value by certain difference.
# I: x (value to add and subtract from), diff (difference)
# O: NULL (if length of x is not 2) or x +/- diff
plusminus = function(x, diff){
  if (length(x) != 2){
    NULL
  } else {
    c(x[1] - diff, x[2] + diff)
  }
}

# Finds the number of pairs that are lower than threshold for sample entropy.
# I: x (RR length data for t-second segment), m (size of group),
#    r (threshold for acceptance), n (size of data)
# O: rate (count / max, proportion of successful groups)
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
  # returns proportion of succesful pairs
  return(count / max)
}

# Converts length type to an index.
# I: length type (length type of observation)
# O: 1, 2, 3 (for "S", "Reg", and "L" respectively)
state_to_ind = function(type){
  return(ifelse(type == "S", 1, ifelse(type == "Reg", 2, 3)))
}

# Determines the index of the transition for the tally of transitions.
# I: prev (length state of previous RR length), 
#    next (length state of current RR length), 
#    counter (tally of transition states)
# O: counter (tally of transition states)
transition_determiner = function(prev, curr, counter){
  prev_ind = state_to_ind(prev)
  curr_ind = state_to_ind(curr)
  ind = 3 * (prev_ind - 1) + curr_ind
  counter[ind] = counter[ind] + 1
  return(counter)
}

# Identifies each observations as AFib, Normal, AFL, or other. Updates 
# tallies.
# I: type (type of observation), counter (tally of types)
# O: counter (tally of types)
type_adder = function(type, counter){
  ind = ifelse(type == "AFIB", 1, ifelse(type == "N", 2, 
                                         ifelse(type == "AFL", 3, 4)))
  counter[ind] = counter[ind] + 1
  return(counter)
}

#### Feature Functions #########################################################

# Finds standardized sample entropy for segment of RR intervals.
# I: x (RR length data for t-second segment), m (size of group),
#    r (threshold for acceptance), n (size of data)
# O: sample entropy (-ln(A / B))
sampEnt = function(x, m = 2, r, n = length(x)){
  # solve for A and B
  A = sampEntCheck(x, m + 1, r, n)
  B = sampEntCheck(x, m, r, n)
  
  # return standardized sample entropy
  return(-log(A / B))
}

# Gets the NEC rate of a segment of intervals.
# I: x (RR lengths), y (RR length differences), int (size of cell),
#    size (number of intervals in segment)
# O: NEC rate (proportion of non-empty cells in segment)
NECRate = function(x, y, int = 0.025, size = length(x)){
  # limits of NEC grid, such that the grid can fit a whole number 
  # of cells of length int
  xl = plusminus(range(x), (int - (diff(range(x)) %% int)) / 2)
  yl = plusminus(range(y), (int - (diff(range(y)) %% int)) / 2)
  
  # grid lines, each separated by int
  xcuts = seq(xl[1], xl[2], by = int); ycuts = seq(yl[1], yl[2], by = int)
  
  # tally of cell counts
  grid = matrix(0, length(ycuts) - 1, length(xcuts) - 1)
  
  # goes through each grouped pair of x and y and 
  # tallies up number of pairs in each cell
  for (p in 1:length(x)){
    # max added in case point is lower bound for x or upper bound for y
    # min added in case point is upper bound for x or lower bound for y
    coords = c(min(max(1, ceiling((x[p] - xl[1]) / int)), dim(grid)[2]), 
               min(max(1, ceiling((yl[2] - y[p]) / int)), dim(grid)[1]))
    
    # update tallies
    grid[coords[2], coords[1]] = grid[coords[2], coords[1]] + 1
  }
  
  # return proportion of non empty cells divided by size (number of intervals)
  return(sum(grid != 0) / size)
}

#### MAIN FUNCTIONS ############################################################

# Function to process each raw atrial fibrillation file of 
# RR intervals for every subject.
# I: number (subject name), 
#    version ("M" for MIT-BIH data, "C" for 2017/CinC data)
# O: Dataset from given file
file_process = function(number, version){
  if (version == "M"){
    number = paste("0", number, sep = "")
    read.csv(paste("Data/MIT-BIH/pre/", number, ".csv", sep = ""), header = T)
  } else {
    read.csv(paste("Data/2017/pre/", number,  ".csv", sep = ""), header = T)
  }
}

# Converts RR intervals into segments, forming new features along the way.
# I: subject (subject name to gather data for),
#    list_data (all RR interval data)
#    source (file path for writing .csv files),
#    chop_left (removes first k indices from whole subject),
#    chop_right (removes last k indices from whole subject),
#    version ("M" or "C"), int (length of segment in seconds),
#    int_min (minimum amount of intervals in segment)
#    r_thresh (multiplier of SD of lengths for sample entropy)
# O: seg_info/NA (segment information for subject)
segment_creator = function(subject, list_data, source, chop_left = 0, 
                           chop_right = 0, version, int, int_min, r_thresh){
  # Gathers data for specific subject
  data = list_data[[subject]]
  
  # AFib/N states to be collected throughout process; added to data frame later
  states = c()
  
  # starting and ending indices for subject under review
  ind_start = 2 + chop_left
  ind_end = dim(data)[1] - chop_right
  
  # tallies up types of transitions (e.g. Reg to Long)
  class_count = rep(0, 9)
  
  # starting index of future segment of t seconds
  seg_start = ind_start - 1
  
  # current starting and ending times for future segment
  start = data$Start[seg_start]
  end = data$End[seg_start]
  
  # gathers lengths and differences
  lengths = c(data$RRLength[seg_start])
  differences = c(data$RRDiff[seg_start])
  
  # used to run sample entropy properly for 2017 data (n >= 5 (int_min))
  if (ind_end - ind_start + 2 < int_min){
    return(NA)
  }
  
  # determines if new data frame is needed
  begin = T
  
  # Goes through each RR interval to create features and segments
  for (i in ind_start:ind_end){
    
    # Creates new end of future segment.
    end = data$End[i]
    
    # Tracks the amount of each type there is in the sample.
    types_count = type_adder(data$State[i], c(0, 0, 0, 0))
    
    # Tracks the lengths and differences of the RR intervals.
    lengths = c(lengths, data$RRLength[i])
    differences = c(differences, data$RRDiff[i])
    
    # Determines which transition occurs based off the RRclass data.
    class_count = transition_determiner(data$RRClass[i - 1], 
                                        data$RRClass[i], class_count)
    
    # For MIT/BIH data: makes new segment when end - start > int
    # For 2017 data: makes new segment when end of truncated interval is reached
    if ((end - start > int & version == "M") | (i == ind_end & version == "C")){
      
      # number of RR intervals inside segment (1 feature)
      intervals = length(lengths)
      
      # Used for MIT-BIH dataset when there is a low number of valid intervals; 
      if (intervals > int_min) {

        # determines if segment is AFIB or not
        state = ifelse(types_count[1] >= sum(types_count[2:4]), "AFIB", "N") 
        states = c(states, state)
        
        # heart rate (1 feature)
        heart_rate = intervals * 60 / (end - start);
        
        # transition proportions (9 features)
        props = round(class_count / (sum(class_count)), digits = 8)
        
        # RR length variance, RR length difference variance, 
        # RR length difference mean (3 features)
        rr_length_var = var(lengths)
        rr_diff_var = var(differences)
        rr_diff_mean = mean(differences)
        
        # sample entropy and NEC rate (2 features)
        samp_ent = sampEnt(x = lengths, m = 2, r = r_thresh * sqrt(rr_length_var))
        necRate = NECRate(x = lengths, y = differences)
        
        # row of information for segment
        row = c(start, end, (end - start), intervals, heart_rate, props, 
                rr_length_var, rr_diff_var, rr_diff_mean, samp_ent, necRate)
        
        # creates new data frame of segment information or adds to data frame
        if (begin == T){
          seg_info = as.data.frame(rbind("1" = row))
          begin = F
        } else {
          seg_info = rbind(seg_info, row)
        }
      } 
      # Resets starting segment time, starting index, types tallies,
      # transition tallies, lengths, and differences for next segment.
      start = data$End[i]
      seg_start = i + 1
      types_count = c(0, 0, 0, 0)
      class_count = rep(0, 9)
      lengths = c()
      differences = c()
    }
  }
  ## DATA WRITING: writes segment info for each subject (MIT-BIH only) and
  ##               returns segment info
  if (begin != T){
    # column names for data frame; binds states information to data frame
    class_names = c("State", "Start", "End", "Time", "Intervals", "HeartRate", 
                    "S-S", "S-Reg", "S-L", "Reg-S", "Reg-Reg", "Reg-L", "L-S", 
                    "L-Reg", "L-L", "RRVar", "dRRVar", "dRRMean", "SampEnt", 
                    "NECRate")
    seg_info = cbind(states, seg_info)
    colnames(seg_info) = class_names
    
    # Writes .csv file with information for particular subject
    if (version == "M"){
      write.csv(seg_info, paste(source, "seg/", subject, ".csv", sep = ""), 
                row.names = FALSE)
    }
    
    # Returns segment information
    return(seg_info)
  } else {
    return(NA)
  }
}

# Combines all subjects' segment information; writes all information to 
# one file. 
# I: list (list of all segment info), source (file path), 
#    version ("M" or "C"), subjects (vector of subjects)
# O: None
combine_to_csv = function(list, source, version, subjects){
  # Used for file naimg purposes.
  vers = ifelse(version == "M", "MIT-BIH", "2017")
  
  # Creates new data frame and adds segment info to data frame.
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
  
  # Writes .csv information.
  write.csv(to_return, 
            paste(source, "seg/all_seg_data_", vers, ".csv", sep = ""), 
            row.names = FALSE)
}

################################################################################