## Jericho Lawson and Georgia Smith
## Summer 2019, 2021
## Adapted from Drew Johnston
## Pre-Processing Program for MIT-BIH Atrial Fibrillation Project (MIT-BIH)

# Using the wfdb package, the code extracts annotations from the MIT-BIH
# Atrial Fibrillation Database (all 23 subjects) and generates basic
# RR-interval information for each of the subjects.

# Note: Make sure your directory has .atr, .dat, .hea, .hea-, and .qrs files
#       for each of the subjects.

#### LIBRARIES ####

import os, wfdb, numpy as np, pandas as pd
from wfdb import processing
from bisect import bisect_left
from STAGE1_MITBIH_functions import outlier_f, run_mean_f, rr_diff_f, rr_class_f 

#### SPECIFICATIONS ####

WEIGHTS = [0.75, 0.25]
RUN_THRESH = [0.85, 1.15]
TYPE_THRESH = "n" # "n" for normal threshold, "q" for quantile threshold
OUT_THRESH = [0.3, 1.5]

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

#### CODE ####

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
    
    # Placeholders for running means, differences in RR intervals, and 
    # transition types.
    rr_mean = [np.nan] * len(rr_int)
    rr_diff = [np.nan] * len(rr_int)
    rr_class = [np.nan] * len(rr_int)
    diffs = 0
    
    # Finds running means, differences in RR intervals, and transition types.
    for obs in range(len(rr_int)):
        # Finds if observation is outlier.
        res = outlier_f(types = TYPE_THRESH, threshes = OUT_THRESH, 
                        length = rr_int[obs], length_data = rr_int)
    
        # Adds three variables to each observation if it is not an outlier.
        # Skips process if it is an outlier.
        if res:
            rr_mean[obs] = run_mean_f(obs, diffs, WEIGHTS, rr_int[obs], 
                                      rr_mean[obs - diffs - 1])
            rr_diff[obs] = rr_diff_f(obs, diffs, rr_int[obs], 
                                     rr_int[obs - diffs - 1])
            rr_class[obs] = rr_class_f(rr_int[obs], rr_mean[obs], RUN_THRESH)
            diffs = 0
        else:
            diffs += 1
    
    # Data formatting.
    data = np.vstack((rr_int, starts, ends, afib, 
                      rr_mean, rr_diff, rr_class)).T
    data = pd.DataFrame(data)
    data.columns = ['RRLength', 'Start', 'End', 'State', 'RRMean', 'RRDiff', 
                    'RRClass']
    data.to_csv(SOURCE_TO + subject + ".csv", index = False)
    
    # Data merging.
    data['Subject'] = [subject] * len(rr_int)
    if subject == '04015':
        combined = data
    else:
        combined = pd.concat([combined, data])
    if subject == '08455':
        combined.to_csv(SOURCE_TO + "all_data_MIT-BIH.csv", index = False)

###############################################################################