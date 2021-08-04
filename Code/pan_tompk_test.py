import os
import matplotlib.pyplot as plt
import numpy as np
import wfdb
from wfdb import processing
import pandas as pd
from bisect import bisect_left
import seaborn as sns 
from pan_tompkin_jericho_comments import pan_tompkin
sns.set()

#os.chdir('C:/Users/mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/mit-bih-atrial-fibrillation-database-1.0.0/files/')
os.chdir('C:/Users/mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Testing/training2017/')

reference = pd.read_csv('REFERENCE.csv',header=None)
reference.columns = ['id','rhythm']
segments = []
zero_rr = []
one_rr = []

for i in range(len(reference.id)):

    subject = reference.id[i]
    rhythm = reference.rhythm[i]
    record = wfdb.rdrecord(subject)

    pan_tompkin(record.p_signal[:,0], record.fs, True, subject)
    # p signal, frame seconds, no plotting, 0 = index of R waves

'''
def spectro(subject):
    signal = wfdb.rdrecord(subject).p_signal[:,0]
    plt.subplot(211)
    plt.xlim(0,len(signal))
    plt.plot(signal)
    plt.subplot(212)
    plt.specgram(signal, Fs=300, cmap='viridis',noverlap=150)
    plt.xlabel('Time (s)')
    plt.ylabel('Frequency (Hz)')
    plt.show()


def peaks_hr(sig, peak_inds, fs, title, figsize=(20, 10), saveto=None):
    "Plot a signal with its peaks and heart rate"
    plt.plot(sig, color='#3979f0', label='Signal')
    plt.plot(peak_inds, sig[peak_inds], 'rx', marker='x', color='#8b0000', label='Peak', markersize=12)
    plt.title(title)
    plt.xlabel('Samples')
    plt.ylabel('ECG (mV)', color='#3979f0')
    #plt.show()

def compare_g_x(name):
    record = wfdb.rdrecord(name)
    qrs_inds = processing.gqrs_detect(sig=record.p_signal[:,0], fs=record.fs) # Compare g to x.
    xqrs_inds = processing.xqrs_detect(sig=record.p_signal[:,0], fs=record.fs, verbose=False, learn=True)
    pt_inds = pan_tompkin(record.p_signal[:,0], record.fs, False)[0]
    plt.subplot(311)
    peaks_hr(sig=record.p_signal, peak_inds=qrs_inds, fs=record.fs, title='GQRS')
    plt.subplot(312)
    peaks_hr(sig=record.p_signal, peak_inds=xqrs_inds, fs=record.fs, title='XQRS')
    plt.subplot(313)
    peaks_hr(sig=record.p_signal, peak_inds=xqrs_inds, fs=record.fs, title='Pan Tompkin')
    plt.tight_layout()
    plt.suptitle(name)
    plt.show()

def drr_plot(subject):
    record = wfdb.rdrecord(subject)
    qrs_inds = processing.gqrs_detect(sig=record.p_signal[:,0], fs=record.fs)
    rr_int = np.array(processing.calc_rr(qrs_inds)/record.fs, dtype='float')
    drr = np.concatenate(([rr_int[0]],processing.calc_rr(rr_int)))
    plt.plot(rr_int, drr, 'o', ms=1)
    plt.title(subject)
    plt.xlabel('RR')
    plt.ylabel('dRR')
    plt.ylim(-.5,.9)
    plt.xlim(.2,1.6)
    plt.show()

GQRS Algo:
zero_rr                                   
['A04735', 'A06103', 'A07183'] (Noise)
No >1.5 seconds: ['A00327','A01133','A01438','A01746','A04735','A04893','A06103','A07183','A07686'] (Other if not Noise from before)
one_rr                                    
[]
No >1.5 seconds: ['A02706', 'A05355'] (Noise)

XQRS Algo:
zero_rr
['A01259','A01445','A01961','A03468','A03549','A03668','A04346','A07528','A08327']
No >1.5 seconds: ['A01133','A01259','A01445','A01961','A02706','A03468','A03549','A03668','A04346','A04893','A05913','A06103','A07528','A07686','A08327']
one_rr
['A01818','A02090','A03044','A03049','A03215','A03260','A04119','A05913','A06103','A06497','A07235']
No >1.5 seconds: ['A01818','A02090','A03044','A03049','A03202','A03215','A03260','A04119','A06497','A07235']


Probably Not Normal, but marked N: (bSQI < 0.7)
0172, 0230, 0316, 0420, 0472, 0818, 0937, 1066, 1148, 1312, 1408, 1445, 1538, 1593, 1713, 1818, 
2023, 2309, 2494, 2513, 2569, 3049, 3066, 3260, 3343, 3412, 3581, 3588, 3623, 3668, 3681, 3918, 
3922, 3974, 4017, 4119, 4338, 4352, 4426, 4804, 4824, 4900, 5324, 5427, 5455, 5555, 5920, 6006, 
6122, 6345, 6457, 6497, 6543, 6697, 6828, 6893, 6895, 7088, 7235, 7416, 7528, 7548, 7721, 7783, 
7933, 8049, 8327

With bSQI < 0.9: 843 Total, 355 are Normal, 267 Other, 67 Afib, 154 Noise
'''