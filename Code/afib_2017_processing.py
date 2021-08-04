import matplotlib.pyplot as plt
import numpy as np
import wfdb
from wfdb import processing
import pandas as pd
from bisect import bisect_left
import seaborn as sns 
from pan_tompkin import pan_tompkin
sns.set()

# Function to calculate running mean.
def running_mean(rr_interval, w1=.75, w2=.25):
    if w1+w2!=1:
        raise ValueError('WEIGHTS MUST ADD TO 1.')
    rmean = [rr_interval[0]]
    for i in range(1,len(rr_interval)):
        rmean.append(w1*rmean[i-1] + w2*rr_interval[i])
    return np.array(rmean).round(decimals=3)
# Funtion to binarize for Transition Matrix, etc.
def binarize(c):
    if c == 'S':
        return 0
    elif c == 'R':
        return 1
    elif c == 'L':
        return 2

# 50 ms grid.
xrange1=np.linspace(.2,1.6,29)
yrange1=np.linspace(-.5,.9,29)
# 25 ms grid.
xrange2=np.linspace(.2,1.6,57)
yrange2=np.linspace(-.5,.9,57)

reference = pd.read_csv('REFERENCE.csv',header=None)
reference.columns = ['id','rhythm']
segments = []
zero_rr = []
one_rr = []
for i in range(len(reference.id)):
    subject = reference.id[i]
    rhythm = reference.rhythm[i]
    record = wfdb.rdrecord(subject)
    qrs_inds = processing.gqrs_detect(sig=record.p_signal[:,0], fs=record.fs)
    # Get bSQI measure.
    xqrs_inds = processing.xqrs_detect(sig=record.p_signal[:,0], fs=record.fs, verbose=False)
    g_count = len(qrs_inds)
    x_count = len(xqrs_inds)
    lower = max(g_count, x_count)
    upper = min(g_count, x_count)
    bSQI = upper/lower
    num_beats = g_count
    seconds = record.sig_len/record.fs
    # Get rr intervals by sample and change to seconds.
    rr_int = np.array(processing.calc_rr(qrs_inds)/record.fs, dtype='float')
    # Don't include data with fewer than 5 heartbeats
    if len(rr_int) <= 5:
        continue
    # Calculate running mean/interval size.
    run_mean = running_mean(rr_int)
    int_size = []
    low_bound = 0.85
    up_bound = 1.15
    for i in range(len(rr_int)):
        if rr_int[i] < low_bound*run_mean[i]:
            int_size.append('S')
        elif rr_int[i] > up_bound*run_mean[i]:
            int_size.append('L')
        else:
            int_size.append('R')
    # Combine and store data.
    drr = np.concatenate(([rr_int[0]],processing.calc_rr(rr_int)))
    ddrr = np.concatenate(([drr[0]],processing.calc_rr(drr)))
    Rvar = drr / run_mean
    data = np.vstack((rr_int, run_mean, drr, ddrr, Rvar, int_size)).T
    data = pd.DataFrame(data)
    data.columns = ['RR','Rmean','dRR','ddRR','Rvar','width']
    data.RR = pd.to_numeric(data.RR) 
    data.Rmean = pd.to_numeric(data.Rmean) 
    data.dRR = pd.to_numeric(data.dRR) 
    data.ddRR = pd.to_numeric(data.ddRR) 
    data.Rvar = pd.to_numeric(data.Rvar)
    ddRR = abs(data.ddRR).mean()
    rvar = data.Rvar.mean()

    #H1 = np.histogram2d(data.RR,data.dRR,bins=[xrange1,yrange1])[0]
    #nec_50ms = np.count_nonzero(H1)
    H2 = np.histogram2d(data.RR,data.dRR,bins=[xrange2,yrange2])[0]
    nec_25ms = np.count_nonzero(H2)
    nec_rate = nec_25ms/num_beats
    # Collect info in one list.
    bpm = 60*num_beats/seconds
    info = [subject, rhythm, num_beats, seconds, bpm, nec_25ms, nec_rate, rvar, ddRR, bSQI]

    binarized = [binarize(c) for c in data.width]
    trans_mat = np.zeros((3,3))
    for (i,j) in zip(binarized,binarized[1:]):
        trans_mat[j,i] += 1
    norm_mat = trans_mat/(len(data.width)-1)
    features = norm_mat.ravel().tolist()

    segments.append(info+features)
full_segs = pd.DataFrame(segments)
full_segs.columns = ['id','rhythm','beats','secs', 'bpm','NEC25ms', 'NECrate','Rvar','ddRR','bSQI','ss','rs','ls','sr','rr','lr','sl','rl','ll']
full_segs.to_csv('2017NEW.csv')

# Divide by afib or not.
af_segs = full_segs[full_segs.rhythm=='A']
normal_segs = full_segs[full_segs.rhythm=='N']
other_segs = full_segs[full_segs.rhythm=='O']
noise_segs = full_segs[full_segs.rhythm=='~']
# Record info
with open('2017NEW.txt','w') as f:
    f.write(full_segs.describe().to_string()+'\n\n')
    f.write('AF:\n'+af_segs.describe().to_string()+'\n\n')
    f.write('Normal:\n'+normal_segs.describe().to_string()+'\n\n')
    f.write('Other:\n'+other_segs.describe().to_string()+'\n\n')
    f.write('Noise:\n'+noise_segs.describe().to_string()+'\n\n')


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