"""
Georgia Smith and Jericho Lawson for NSF Grant #1659288

This is a Python File that will include functions that can be accessed and called by R
to extract data from Physionet using the python WFDB package which is not available from CRAN

Dependencies can be installed with one of the following commands [Given that python2 is already installed]
  - pip install numpy matplotlib sklearn wfdb scipy pandas tqdm
  - conda install numpy matplotlib scikit-learn wfdb scipy pandas tqdm
"""

import numpy as np
from matplotlib import pyplot as plt
from sklearn import preprocessing
import wfdb
import copy as cp
import scipy.signal as signal
from sklearn import preprocessing
from tqdm import tqdm
import os
import re
import pandas as pd
import pickle
import csv
import requests, zipfile, io
import warnings

def download_mit_bih(download_path = "../mit-bih-raw/"):
    if os.path.exists(download_path+"/files"):
        raise SystemError('MIT-BIH data already downloaded at ' + download_path)
        return
    
    mit_bih_raw = requests.get("https://physionet.org/static/published-projects/afdb/mit-bih-atrial-fibrillation-database-1.0.0.zip")
    mit_bih_zipfile = zipfile.ZipFile(io.BytesIO(mit_bih_raw.content))
    
    if os.path.isdir(download_path):
        mit_bih_zipfile.extractall(download_path)
    else: 
        os.mkdir(download_path)
        mit_bih_zipfile.extractall(download_path)



def extract_mit_bih(mitbih_path = '../mit-bih-raw/files/', save_path = "../mit_bih_extracted/", reload_flag = True):
    rlist = []
    records = mitbih_path + "RECORDS"
    if not os.path.exists(records):
        sys.exit("\n\tERROR: path " + records + " does not exist!\n\tRun function download_mit_bih() or check file location for the RECORDS file")
        return
    
    if not os.path.isdir(save_path):
        os.mkdir(save_path)
        os.mkdir(save_path + "dataframes/")
        os.mkdir(save_path + "extracted_csv/")
    elif not os.path.isdir(save_path + "dataframes/"):
        os.mkdir(save_path + "dataframes/")
        if not os.path.isdir(save_path + "extracted_csv/"):
            os.mkdir(save_path + "extracted_csv/")
    elif not os.path.isdir(save_path + "extracted_csv/"):
        os.mkdir(save_path + "extracted_csv/")

    with open(records) as rfile:
        for record in rfile:
            record = record[0:len(record)-1]
            rlist.append(record)
            
    ###### Step 1: Initialize all Arrays
    samples = [] # will house the samples of all subjects
    good_list = [] # will list the names of the subjects we successfully extracted
    bad_list = [] # will house the names of the subjects we failed to extract
    qrs = [] # will house the indices of R-Peaks for all subjects
    atr_label = [] # will house the labels for each rhythm annotation for all subjects
    atr_locs = [] # will house the locations corresponding to the rhythm annotation labels


    ###### Step 2: Extract Information
    for x in tqdm(rlist):
        try: # A try statement will run the except statement if for some reason the try commands fail
             # In this case I use the try statement because one of the subjects has no signal data causing failure
             # I then use bad_list and good_list so that all of the indices in rlist match with the arrays we initialized in Step 1, above
            ######################################################
            samp = wfdb.rdsamp(mitbih_path+x) # wfdb.rdsamp(file_location) will read the signal & header data and return a 2 value array
                # samp[0] - the signal data is the raw reading from the ecg. Each value is a sample taken.
                # samp[1] - the header data includes things about the signal data such as:
                  # samples per section, denoted 'fs'
                  # number of signals, denoted 'n_sig'
            ######################################################
            samples.append(samp) #add it to our array for all subject

            qrs_tmp = wfdb.rdann(mitbih_path+x, extension='qrs') #extract the QRS Info
            qrs_locs = np.array(qrs_tmp.sample, dtype='int') #Get just the loccation of R-Peaks from the QRS Info
            qrs.append(qrs_locs) # Add to our array for all subjects

            atr = wfdb.rdann(mitbih_path+x,extension='atr') #extract the atr info which stores the rhythm type(s) over the whole signal
            atr_label.append(atr.aux_note) # aux_note stores the type of rhythm - main two are '(N' for normal and '(AFIB' for AFIB
            atr_locs.append(np.append(atr.sample, len(samp[0]))) #I add the length of the whole sample to the end for better visualization later

            good_list.append(x) # when all extraction is successful append the record name to good_list
        except Exception as exep:
            print(exep) # Alert the user of an exception
            bad_list.append(x) # add to the bad list
    
    subject_dataframes = [] # Initialize the subject_dataframes - will hold all of our subject dataframes

    atr_dics = [] #Initialize the array that will hold the dictionary for each subject

    for idxs,lab in enumerate(atr_label):
        atr_dic = {} #Initialize dictionary for each subject
        for idx,x in enumerate(lab):
            if x not in atr_dic.keys():
                atr_dic[x] = [] #Add dictionary key if does not exist
            atr_dic[x].append([atr_locs[idxs][idx], atr_locs[idxs][idx+1]]) #Insert range for each rhythm
        atr_dics.append(atr_dic) #Add to dictionary array
    
    
    for s, _ in enumerate(tqdm(good_list)): # Iterate through all of the subjects that we have complete data of 
        subj = pd.DataFrame( # The below statements initialize our datafram. The first to columns will be our given signals, and the rest we initialize to 0
            data = np.transpose(np.array([ # First we give our data, for pandas they want the data by row instead of by column, so we use transpose to get the proper format
                                                   [x[0] for x in samples[s][0]],
                                                   [x[1] for x in samples[s][0]],
                                                   np.zeros(len(samples[s][0])), # np.zeros makes an array of zeros with the given lenth
                                                   np.zeros(len(samples[s][0])), 
                                                   np.zeros(len(samples[s][0])), 
                                                   np.zeros(len(samples[s][0])), 
                                            ])
                               ),
            columns = ['Signal 1', 'Signal 2', 'R-Peak', 'Normal', 'AFIB', 'Other'] # Here we name our columns to match the dataframe we outlined above
        )
        norm = [] # Initialize the norm array which will list every index the person is in a normal rhythm
        if '(N' in atr_dics[s].keys():
            for x in atr_dics[s]['(N']: # Then we iterate through our ranges we extracted above
                norm = norm + list(range(x[0], x[1])) # And add all values in the range to our norm array
        af = [] # Then we do the same steps above for AFIB rhythms
        if '(AFIB' in atr_dics[s].keys():
            for x in atr_dics[s]['(AFIB']:
                af = af + list(range(x[0], x[1]))
        subj['R-Peak']= subj.index.isin(qrs[s]) # the isin() function of a DataFram index will return true if the index is in that list and false if it is not
                                                # then, we can initialize our dataFrame with correct values based on that
        subj['Normal']= subj.index.isin(norm)
        subj['AFIB'] = subj.index.isin(af)
        subj['Other'] = ~subj.index.isin(np.append(norm, af)) # Because we are classifying AFIB specifically we define other as any rhythm not in the norm or AFIB list

        subject_dataframes.append(subj) # Add the dataframe we built to our to array that holds all of our subjects' dataframes
        
    for idx, x in enumerate(tqdm(good_list)): 
        subject_dataframes[idx].to_csv(save_path + 'dataframes/'+x+'.csv') # Pandas DataFrames have a built in to_csv() function which whill save it at the passed path

    np.savetxt(save_path + 'dataframes/' + "subject_list.csv", good_list, delimiter=",",  fmt='%s') 
       # We'll load the complete list of subjects as well so that we can easily recreate the file names

    np.savetxt(save_path + "extracted_csv/" + "subject_list.csv", good_list, delimiter=",",  fmt='%s') #Save the names in the folder 
    for idx, x in enumerate(tqdm(good_list)): # Iterate through our subjects
        if not os.path.exists(save_path + "extracted_csv/"+x+"_signals.csv") or reload_flag:
            np.savetxt(save_path + "extracted_csv/"+x+"_signals.csv", np.array(samples[idx][0]), delimiter=",") # numPy has a savetxt() function which by setting the delimiter as ',' we can 
                                                                                                # simulate a to_csv() function 
        if not os.path.exists(save_path + "extracted_csv/"+x+"_rpeaks.csv") or reload_flag:
                np.savetxt(save_path + "extracted_csv/"+x+"_rpeaks.csv", np.array(qrs[idx]), delimiter=",")      
        if not os.path.exists(save_path + "extracted_csv/"+x+"_headers.pkl") or reload_flag:
            with open(save_path + "extracted_csv/"+x+"_headers.pkl", 'wb') as picklefile: # nomPy has no way to save a dictionary as a CSV so we use the pickle package
                                        # First we open up the file we would like to write to
                pickle.dump(samples[idx][1], picklefile)
        if not os.path.exists(save_path + "extracted_csv/"+x+"_rhythms.pkl") or reload_flag:
            with open(save_path + "extracted_csv/"+x+"_rhythms.pkl", 'wb') as picklefile:
                pickle.dump(atr_dics[idx], picklefile)

