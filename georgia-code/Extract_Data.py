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



def extract_mit_bih(mitbih_path = '../mit-bih-raw/files/'):
    rlist = []
    records = mitbih_path + "RECORDS"
    if not os.path.exists(records):
        sys.exit("\n\tERROR: path " + records + " does not exist!\n\tRun function download_mit_bih() or check file location for the RECORDS file")
        return

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
    
    return(list(good_list), list(bad_list), list(qrs), list(atr_label), list(atr_locs))
          