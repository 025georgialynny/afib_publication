## Jericho Lawson
## Summer 2019, 2021
## AFib Features

# Used to create features for AFib classification, using the STAGE2functions.R file.

#### SPECIFICATIONS ####
INTERVAL = 30
VERSION = "C" # M for MIT-BIH, C for 2017 CinC dataset
DIR = "C:/Users/mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Summer 2021 Testing/" # directory with RR data files

#### CODE ####
# code with functions to run processes below
setwd(DIR)
source("Code/STAGE2functions.R")

# Retrieves all R-R interval information for all subjects.
if (VERSION == "M"){ # "Data/MIT-BIH" for MIT-BIH, "Data/2017" for 2017 CinC dataset
  subpath = "Data/MIT-BIH/"
  subjects = as.character(read.csv(paste(subpath, "ref_MIT-BIH.csv", sep = ""), header = F, fileEncoding = "UTF-8-BOM")$V1)
} else {
  subpath = "Data/2017/"
  subjects = read.csv(paste(subpath, "ref_2017.csv", sep = ""), header = F)$V1
}


data = lapply(subjects, file_process, version = VERSION)
data = lapply(data, omitter)

names(data) = subjects

data2 = lapply(subjects, segment_creator, list_data = data, source = paste(DIR, subpath, sep = ""), version = VERSION, int = INTERVAL, int_min = 5)
combine_to_csv(data2, source = paste(DIR, subpath, sep = ""), version = VERSION, subjects = subjects)

###################################################################################################