## Jericho Lawson
## Summer 2019, 2021
## AFib Features

# Used to create features for AFib classification, 
# using the STAGE2functions.R file.

#### SPECIFICATIONS ####################################################################################################

INTERVAL = 30 # length of segments in seconds (MIT-BIH data only)
VERSION = "C" # M for MIT-BIH, C for 2017 CinC dataset
COEF_SE = 2 # multiplied by standard deviation for sample entropy
WRITE_SEP = F # writes individual and combined .csv files if true; writes combined file otherwise
# directory with RR data files
DIR = "C:/Users/mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Summer 2021 Testing/" 

#### CODE ##############################################################################################################

# Sets directory; gathers file with functions to run program.
setwd(DIR)
source("Code/STAGE2functions.R")

# Retrieves all R-R interval information for all subjects.
if (VERSION == "M"){ # "Data/MIT-BIH" for MIT-BIH, "Data/2017" for 2017 CinC dataset
  subpath = "Data/MIT-BIH/"
  subjects = as.character(read.csv(paste(subpath, "ref_MIT-BIH.csv", sep = ""), 
                                   header = F, fileEncoding = "UTF-8-BOM")$V1)
} else {
  subpath = "Data/2017/"
  subjects = read.csv(paste(subpath, "ref_2017.csv", sep = ""), header = F)$V1
}

# Processes .csv files that contain RR interval data; omits any intervals 
# that have NA data due to being too short or long
data = lapply(subjects, file_process, version = VERSION)
data = lapply(data, na.omit)

# Adds subjects' names to list of data
names(data) = subjects

# Uses apply function to create dataset of segments with various features; 
# combines all data together to make one .csv file
data2 = lapply(subjects, segment_creator, list_data = data, 
               source = paste(DIR, subpath, sep = ""), 
               version = VERSION, int = INTERVAL, int_min = 5, r_thresh = COEF_SE, write_sep = WRITE_SEP)
combine_to_csv(data2, source = paste(DIR, subpath, sep = ""), 
               version = VERSION, subjects = subjects, type = INTERVAL)

############################################################################################################################