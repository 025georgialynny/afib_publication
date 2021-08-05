## Jericho Lawson
## Summer 2019, 2021
## Merging and Organization

# Adds calculated features to the data gathered from the MIT-BIH/2017 subjects.
# Merges all data into one file.

## Set directory.
DIR = "C:/Users/mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Summer 2021 Testing/"
setwd(DIR)

## Pick the correct version for testing. ##
VERSION = "C" # M for MIT-BIH, C for 2017 CinC dataset

# Array of subjects.
if (VERSION == "M"){
  subjects = as.character(read.csv("Data/MIT-BIH/ref_MIT-BIH.csv", header = F, fileEncoding = "UTF-8-BOM")$V1)
  setwd(paste(DIR, "Data/MIT-BIH/pre/", sep = ""))
} else {
  subjects = read.csv("Data/2017/ref_2017.csv", header = F)$V1
  setwd(paste(DIR, "Data/2017/pre/", sep = ""))
}

# Reads in data, adds subject column, and merged data together.
for (spot in 1:length(subjects)){
  sub = ifelse(VERSION == "M", paste("0", subjects[spot], sep = ""), subjects[spot])
  new = read.csv(paste(sub, ".csv", sep = ""), header = T)
  new$Subject = subjects[spot]
  if (spot == 1){
    merged = new
  } else {
    merged = rbind(merged, new)
  }
}


# Write file.
if (VERSION == "M"){
  write.csv(merged, "all_data_MIT-BIH.csv", row.names = FALSE)
} else {
  write.csv(merged, "all_data_2017.csv", row.names = FALSE)
}

###################################################################################################