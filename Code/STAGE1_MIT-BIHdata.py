## Jericho Lawson
## Summer 2019, 2021
## Adapted from Drew Johnston
## Pre-Processing Program for MIT-BIH Atrial Fibrillation Project

# Using the wfdb package, the code extracts annotations from the MIT-BIH
# Atrial Fibrillation Database (all 23 subjects) and generates basic
# RR-interval information for each of the subjects.

# Note: Make sure your directory has .atr, .dat, .hea, .hea-, and .qrs files
#       for each of the subjects.

#### LIBRARIES ####
import os, wfdb, numpy as np, pandas as pd
from wfdb import processing
from bisect import bisect_left

# Directory where all files are located
SOURCE_FROM = "C:/Users/mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/mit-bih-atrial-fibrillation-database-1.0.0/files/"
os.chdir(SOURCE_FROM)

# Directory where new .csv files will be located
SOURCE_TO = "C:/Users/mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Summer 2021 Testing/Data/MIT-BIH/pre/"
# Note: 00735 and 03665 not included due to incompleteness.
FILES = ['04015', '04043', '04048', '04126', '04746', '04908', '04936',
         '05091', '05121', '05261', '06426', '06453', '06995', '07162',
         '07859', '07879', '07910', '08215', '08219', '08378', '08405',
         '08434', '08455']

for subject in FILES:
    # Extracts QRS and ANN information.
    qrs = wfdb.rdann(subject, 'qrs') # qrs.sample has heartbeats
    ann = wfdb.rdann(subject, 'atr') # ann.aux_note has types, ann.sample has time of changes in type
    
    # Identifies type during each heartbeat. 
    # Note: heart rate ends with the previous type.
    afib = []
    for x in qrs.sample[1:] - 1:
        i = bisect_left(ann.sample, x)
        if i == 0:
            afib.append("N")
        else:
            afib.append(ann.aux_note[i - 1][1:])
            
    # Identifies starting and ending times, as well as RR interval lengths.
    starts = []; ends = []
    for index in range(len(qrs.sample) - 1):
        starts.append(qrs.sample[index] / qrs.fs);
    ends = starts[1:]; ends.append(qrs.sample[len(qrs.sample) - 1] / qrs.fs)
    rr_int = processing.calc_rr(qrs.sample)/ qrs.fs # divides by 250 (fs)
    
    # Data formatting.
    data = np.vstack((rr_int, starts, ends, afib)).T
    data = pd.DataFrame(data)
    data.columns = ['RRLength', 'Start', 'End', 'State']
    data.to_csv(SOURCE_TO + subject + ".csv", index = False)
    
###############################################################################