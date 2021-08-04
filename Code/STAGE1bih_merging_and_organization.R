## Jericho Lawson
## Summer 2019, 2021
## Merging and Organization

# Adds calculated features to the data gathered from the MIT-BIH/2017 subjects.
# Merges all data into one file.

## Pick the correct version for testing. ##
VERSION = "M"
#VERSION = "C"
#VERSION = "multi_2017"

# Array of subjects.
if (VERSION == "M"){
  subjects = c('04015', '04043', '04048', '04126', '04746', '04908', '04936', '05091', '05121', 
               '05261', '06426', '06453', '06995', '07162', '07859', '07879', '07910', '08215', 
               '08219', '08378', '08405', '08434', '08455')
}else{
  subjects = read.csv("C:/Users/mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Testing/training2017/REFERENCE_2.csv", header = F)$V2[-1]
}

# Reads in data, adds subject column, and merged data together.
for (spot in 1:length(subjects)){
  if (VERSION == "M"){
    new = read.csv(paste("C:/Users/mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Data/", subjects[spot], ".csv", sep = ""), header = T)
  }else if (VERSION == "C"){
    new = read.csv(paste("C:/Users/mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Data/Training_Set/", subjects[spot], ".csv", sep = ""), header = T)
  }else{
    new = read.csv(paste("C:/Users/mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Data/Training_Set/", subjects[spot], "_multi.csv", sep = ""), header = T)
  }
  new$Subject = subjects[spot]
  if (spot == 1){
    merged = new
  }else{
    merged = rbind(merged, new)
  }
}

# Write file.
if (VERSION == "M"){
  write.csv(merged, "C:/Users/mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Data/all_data.csv")
}else if (VERSION == "C"){
  write.csv(merged, "C:/Users/mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Data/all_data_2017_train.csv")
}else{
  write.csv(merged, "C:/Users/mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Data/all_data_2017_train_multi.csv")
}
write.csv(merged, "C:/Users/mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Data/all_data.csv")

###################################################################################################