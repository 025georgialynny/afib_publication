## Jericho Lawson
## Summer 2019
## Atrial Fibrillation - Testing: Text to CSV for 2017 Data
# NOT NEEDED ANYMORE

#### Library for cleaning testing 2017 data. ####
library(tidyr)

#### Data Choice? ####
VERSION = "single_2017"
#VERSION = "multi_2017"

#### Collect testing 2017 state data. ####
if (VERSION == "single_2017"){
  states = read.csv("C:/Users/mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Testing/sample2017/validation/REFERENCE.csv", header = F)
  states$V2 = ifelse(states$V2 != "A", "N", "AFIB") # use for binary classification
}else if (VERSION == "multi_2017"){
  states = read.csv("C:/Users/mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Testing/training2017/REFERENCE.csv", header = F)
}

#### Collect data. ####

# Collect testing 2017 RRLength data, merge with state data, and write new csv 
# file for each subject.
if (VERSION == "single_2017"){
  for (num in 1:length(states$V2)){
    if (states$V1[num] != "A04735"){
      new_data = read.delim2(paste("C:/Users/mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Testing/sample2017/", 
                                  states$V1[num], ".txt", sep = ""), header = F)
      new_data = separate(new_data, col = V1, sep = ",", into = c("Start", "End", "RRLength"))
      new_data$State = paste("(", states$V2[num], sep = "")
      write.csv(new_data, paste("C:/Users/mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Data/Testing_Set/", 
                                states$V1[num], ".csv", sep = ""))
    }
  }
}else{ #### For training 2017 data.
  aa = c()
  #### Collect testing 2017 RRLength data, merge with state data, and write new csv file for 
  #### each subject.
  for (num in 1:length(states$V2)){
    if (file.info(paste("C:/Users/mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Testing/training2017/", 
                        states$V1[num], ".txt", sep = ""))$size != 0){
      aa = c(aa, paste(states$V1[num]))
      new_data = read.delim2(paste("C:/Users/mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Testing/training2017/", 
                                   states$V1[num], ".txt", sep = ""), header = F)
      new_data = separate(new_data, col = V1, sep = ",", into = c("Start", "End", "RRLength"))
      new_data$State = paste("(", states$V2[num], sep = "")
      write.csv(new_data, paste("C:/Users/mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Data/Training_Set/", 
                                states$V1[num], "_multi.csv", sep = ""))
    }
  }
  write.csv(aa, paste("C:/Users/mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Testing/training2017/REFERENCE_2.csv", sep = ""))
}
###################################################################################################