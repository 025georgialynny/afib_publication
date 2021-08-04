#Vignesh Kalidas 
#PhD - Computer Engineering, 
#Dept of Electrical Engineering, 
#University of Texas at Dallas, Texas, USA
# GNU General public License
# Complete implementation of Pan-Tompkins algorithm

# Translated to Python 3.6 and adapted by Drew Johnston

## Inputs
# ecg : raw ecg vector signal 1d signal
# fs : sampling frequency e.g. 200Hz, 400Hz and etc
# gr : flag to plot or not plot 

## Outputs
# qrs_amp_raw : amplitude of R waves amplitudes
# qrs_i_raw : index of R waves
# delay : number of samples which the signal is delayed due to the filtering

## Method :

## PreProcessing
# 1) Signal is preprocessed , if the sampling frequency is higher then it is downsampled
# and if it is lower upsampled to make the sampling frequency 200 Hz
# with the same filtering setups introduced in Pan
# tompkins paper (a combination of low pass and high pass filter 5-15 Hz)
# to get rid of the baseline wander and muscle noise. 

# 2) The filtered signal
# is derivated using a derivating filter to high light the QRS complex.

# 3) Signal is squared.

# 4) Signal is averaged with a moving window to get rid
# of noise (0.150 seconds length).

# 5) depending on the sampling frequency of your signal the filtering
# options are changed to best match the characteristics of your ecg signal

# 6) Unlike the other implementations in this implementation the desicion
# rule of the Pan tompkins is implemented completely.

## Decision Rule 
# At this point in the algorithm, the preceding stages have produced a roughly pulse-shaped
# waveform at the output of the MWI . The determination as to whether this pulse
# corresponds to a QRS complex (as opposed to a high-sloped T-wave or a noise artefact) is
# performed with an adaptive thresholding operation and other decision
# rules outlined below

# a) FIDUCIAL MARK - The waveform is first processed to produce a set of weighted unit
# samples at the location of the MWI maxima. This is done in order to localize the QRS
# complex to a single instant of time. The w[k] weighting is the maxima value.

# b) THRESHOLDING - When analyzing the amplitude of the MWI output, the algorithm uses
# two threshold values (THR_SIG and THR_NOISE, appropriately initialized during a brief
# 2 second training phase) that continuously adapt to changing ECG signal quality. The
# first pass through y[n] uses these thresholds to classify the each non-zero sample
# (CURRENTPEAK) as either signal or noise:
# If CURRENTPEAK > THR_SIG, that location is identified as a �QRS complex
# candidate� and the signal level (SIG_LEV) is updated:
# SIG _ LEV = 0.125 �CURRENTPEAK + 0.875� SIG _ LEV

# If THR_NOISE < CURRENTPEAK < THR_SIG, then that location is identified as a
# �noise peak� and the noise level (NOISE_LEV) is updated:
# NOISE _ LEV = 0.125�CURRENTPEAK + 0.875� NOISE _ LEV
# Based on new estimates of the signal and noise levels (SIG_LEV and NOISE_LEV,
# respectively) at that point in the ECG, the thresholds are adjusted as follows:
# THR _ SIG = NOISE _ LEV + 0.25 � (SIG _ LEV ? NOISE _ LEV )
# THR _ NOISE = 0.5� (THR _ SIG)
# These adjustments lower the threshold gradually in signal segments that are deemed to
# be of poorer quality.


# c) SEARCHBACK FOR MISSED QRS COMPLEXES - In the thresholding step above, if
# CURRENTPEAK < THR_SIG, the peak is deemed not to have resulted from a QRS
# complex. If however, an unreasonably long period has expired without an abovethreshold
# peak, the algorithm will assume a QRS has been missed and perform a
# searchback. This limits the number of false negatives. The minimum time used to trigger
# a searchback is 1.66 times the current R peak to R peak time period (called the RR
# interval). This value has a physiological origin - the time value between adjacent
# heartbeats cannot change more quickly than this. The missed QRS complex is assumed
# to occur at the location of the highest peak in the interval that lies between THR_SIG and
# THR_NOISE. In this algorithm, two average RR intervals are stored,the first RR interval is 
# calculated as an average of the last eight QRS locations in order to adapt to changing heart 
# rate and the second RR interval mean is the mean 
# of the most regular RR intervals . The threshold is lowered if the heart rate is not regular 
# to improve detection.

# d) ELIMINATION OF MULTIPLE DETECTIONS WITHIN REFRACTORY PERIOD - It is
# impossible for a legitimate QRS complex to occur if it lies within 200ms after a previously
# detected one. This constraint is a physiological one � due to the refractory period during
# which ventricular depolarization cannot occur despite a stimulus[1]. As QRS complex
# candidates are generated, the algorithm eliminates such physically impossible events,
# thereby reducing false positives.

# e) T WAVE DISCRIMINATION - Finally, if a QRS candidate occurs after the 200ms
# refractory period but within 360ms of the previous QRS, the algorithm determines
# whether this is a genuine QRS complex of the next heartbeat or an abnormally prominent
# T wave. This decision is based on the mean slope of the waveform at that position. A slope of
# less than one half that of the previous QRS complex is consistent with the slower
# changing behaviour of a T wave � otherwise, it becomes a QRS detection.
# Extra concept : beside the points mentioned in the paper, this code also
# checks if the occured peak which is less than 360 msec latency has also a
# latency less than 0,5*mean_RR if yes this is counted as noise

# f) In the final stage , the output of R waves detected in smoothed signal is analyzed and double
# checked with the help of the output of the bandpass signal to improve
# detection and find the original index of the real R waves on the raw ecg
# signal

import numpy as np 
import matplotlib.pyplot as plt 
from scipy.signal import lfilter, convolve, butter, filtfilt, find_peaks
import wfdb
from wfdb import processing
import pandas as pd

def pan_tompkin(ecg, fs=300, gr=True, subject_name='Unknown'):
    if ecg.ndim != 1:
        raise ValueError('ecg must be a row or column vector')

    # Initialize
    qrs_c =[] #amplitude of R
    qrs_i =[] #index
    SIG_LEV = 0
    nois_c =[]
    nois_i =[]
    delay = 0
    skip = 0 # becomes one when a T wave is detected
    not_nois = 0 # it is not noise when not_nois = 1
    selected_RR =[] # Selected RR intervals
    m_selected_RR = 0
    mean_RR = 0
    qrs_i_raw =[]
    qrs_amp_raw=[]
    ser_back = 0
    test_m = 0
    SIGL_buf = []
    NOISL_buf = []
    THRS_buf = []
    SIGL_buf1 = []
    NOISL_buf1 = []
    THRS_buf1 = []


    ## Plot differently based on filtering settings
    if gr:
        ax1 = plt.subplot(321)
        plt.tight_layout()
        plt.title("Raw ECG Signal")
        plt.plot(ecg)

    ## Noise cancelation(Filtering) # Filters (Filter in between 5-15 Hz)
    if fs == 200:
        ## Low Pass Filter  H(z) = ((1 - z^(-6))^2)/(1 - z^(-1))^2
        b = np.array([1, 0, 0, 0, 0, 0, -2, 0, 0, 0, 0, 0, 1])
        a = np.array([1, -2, 1])
        h_l = lfilter(b, a, [1.]+[0.]*12, axis=0)
        ecg_l = np.convolve(ecg, h_l)
        ecg_l = ecg_l/max(abs(ecg_l))
        delay = 6 #based on the paper
        if gr:
            ax2 = plt.subplot(322, sharey=ax1)
            plt.title('Low Pass Filtered')
            plt.plot(ecg_l)
        ## High Pass filter H(z) = (-1+32z^(-16)+z^(-32))/(1+z^(-1))
        b = [-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 32, -32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
        a = [1, -1]
        h_h = lfilter(b, a, [1.]+[0.]*32, axis=0)
        ecg_h = np.convolve(ecg_l, h_h)
        ecg_h = ecg_h/max(abs(ecg_h))
        delay = delay + 16 # 16 samples for highpass filtering
        if gr:
            ax3 = plt.subplot(323)
            plt.title('High Pass Filtered')
            plt.plot(ecg_h)

    else:
        if gr:
            ax2 = plt.subplot(322)
            plt.title("Derivative")
            plt.plot(np.gradient(ecg))
        ## bandpass filter for Noise cancelation of other sampling frequencies(Filtering)
        f1 = 5 #cuttoff low frequency to get rid of baseline wander
        f2 = 15 #cuttoff frequency to discard high frequency noise
        Wn = np.array([f1, f2])*2/fs # cutt off based on fs
        N = 3 # order of 3 less processing
        a,b = butter(N,Wn,btype='bandpass') #bandpass filtering
        ecg_h = filtfilt(a,b,ecg)
        ecg_h = ecg_h/max(abs(ecg_h))
        if gr:
            ax3 = plt.subplot(323)
            plt.title("Band Pass Filtered")
            plt.plot(ecg_h)
    ## derivative filter H(z) = (1/8T)(-z^(-2) - 2z^(-1) + 2z + z^(2))
    h_d = np.array([-1., -2, 0, 2, 1])*(1/8)#1/8*fs
    ecg_d = np.convolve(ecg_h ,h_d)
    ecg_d = ecg_d/max(ecg_d)
    delay = delay + 2 # delay of derivative filter 2 samples
    if gr:
        ax4 = plt.subplot(324)
        plt.title("Filtered with the derivative filter")
        plt.plot(ecg_d)
    
    ## Squaring nonlinearly enhance the dominant peaks
    ecg_s = ecg_d**2
    if gr:
        ax5 = plt.subplot(325)
        plt.title("Squared")
        plt.plot(ecg_s)

    ## Moving average Y(nt) = (1/N)[x(nT-(N - 1)T)+ x(nT - (N - 2)T)+...+x(nT)]
    ecg_m = np.convolve(ecg_s ,np.ones(round(0.150*fs))/round(0.150*fs))
    delay = delay + 15
    if gr:
        ax6 = plt.subplot(326)
        plt.title("Averaged with 30 samples")
        # length,Black noise,Green Adaptive Threshold,RED Sig Level,
        # Red circles QRS adaptive threshold
        plt.plot(ecg_m)

    
    ## Fiducial Mark 
    # Note : a minimum distance of 40 samples is considered between each R wave
    # since in physiological point of view no RR wave can occur in less than
    # 200 msec distance
    #[pks,locs] = findpeaks(ecg_m,'MINPEAKDISTANCE',round(0.2*fs))
    locs, properties = find_peaks(ecg_m, distance=round(0.2*fs), height=.001)
    peaks = properties['peak_heights']

    ## initialize the training phase (2 seconds of the signal) to determine the THR_SIG and THR_NOISE
    if len(ecg_m) <= 2*fs: # Length of the signal too small
        THR_SIG = max(ecg_m)*1/4 # 0.25 of the max amplitude 
        THR_NOISE = np.mean(ecg_m)*1/2 # 0.5 of the mean signal is considered to be noise
    else:
        THR_SIG = max(ecg_m[:2*fs])*1/4 # 0.25 of the max amplitude 
        THR_NOISE = np.mean(ecg_m[:2*fs])*1/2 # 0.5 of the mean signal is considered to be noise
    SIG_LEV = THR_SIG
    NOISE_LEV = THR_NOISE


    ## Initialize bandpath filter threshold(2 seconds of the bandpass signal)
    if len(ecg_h) <= 2*fs: # If length of the signal is too small
        THR_SIG1 = max(ecg_h)*1/4 # 0.25 of the max amplitude 
        THR_NOISE1 = np.mean(ecg_h)*1/2 #
    else:
        THR_SIG1 = max(ecg_h[:2*fs])*1/4 # 0.25 of the max amplitude 
        THR_NOISE1 = np.mean(ecg_h[:2*fs])*1/2 #
    SIG_LEV1 = THR_SIG1 # Signal level in Bandpassed filter
    NOISE_LEV1 = THR_NOISE1 # Noise level in Bandpassed filter
    
    ## Thresholding and online decision rule
    for i in range(len(peaks)):
    ## locate the corresponding peak in the filtered signal 
        if locs[i]-round(0.150*fs) >= 0 and locs[i] <= len(ecg_h):
            y_i = max(ecg_h[locs[i]-round(0.150*fs):locs[i]])
            x_i = np.argmax(ecg_h[locs[i]-round(0.150*fs):locs[i]])
        else:
            if i == 0:
                y_i = max(ecg_h[:locs[i]])
                x_i = np.argmax(ecg_h[:locs[i]])
                ser_back = 1
            elif locs[i] >= len(ecg_h):
                y_i = max(ecg_h[locs[i]-round(0.150*fs):])
                x_i = np.argmax(ecg_h[locs[i]-round(0.150*fs):])


    ## update the heart_rate (Two heart rate means one the most recent and the other selected)
        if len(qrs_c) >= 9:
            diffRR = processing.calc_rr(qrs_i[-9:]) #calculate RR interval
            mean_RR = np.mean(diffRR) # calculate the mean of 8 previous R waves interval
            comp = qrs_i[-1]-qrs_i[-2] #latest RR
            if comp <= 0.92*mean_RR or comp >= 1.16*mean_RR:
                # lower down thresholds to detect better in MVI
                    THR_SIG = 0.5*THR_SIG
                    #THR_NOISE = 0.5*(THR_SIG)  
                # lower down thresholds to detect better in Bandpass filtered 
                    THR_SIG1 = 0.5*THR_SIG1
                    #THR_NOISE1 = 0.5*(THR_SIG1) 
            else:
                m_selected_RR = mean_RR #the latest regular beats mean

            ## calculate the mean of the last 8 R waves to make sure that QRS is not
            # missing(If no R detected , trigger a search back) 1.66*mean

        if m_selected_RR:
            test_m = m_selected_RR #if the regular RR availabe use it   
        elif mean_RR and m_selected_RR == 0:
            test_m = mean_RR   
        else:
            test_m = 0

        if test_m:
            if (locs[i] - qrs_i[-1]) >= round(1.66*test_m):# it shows a QRS is missed 
                pks_temp = max(ecg_m[qrs_i[-1] + round(0.2*fs):locs[i] - round(0.2*fs)]) # search back and locate the max in this interval
                locs_temp = np.argmax(ecg_m[qrs_i[-1] + round(0.2*fs):locs[i] - round(0.2*fs)])
                locs_temp = qrs_i[-1] + round(0.2*fs) + locs_temp - 1 #location 

                if pks_temp > THR_NOISE:
                    qrs_c.append(pks_temp)
                    qrs_i.append(locs_temp)
                    # find the location in filtered sig
                    if locs_temp <= len(ecg_h):
                        y_i_t = max(ecg_h[locs_temp-round(0.150*fs):locs_temp])
                        x_i_t = np.argmax(ecg_h[locs_temp-round(0.150*fs):locs_temp])
                    else:
                        y_i_t = max(ecg_h[locs_temp-round(0.150*fs):])
                        x_i_t = np.argmax(ecg_h[locs_temp-round(0.150*fs):])
                    # take care of bandpass signal threshold
                    if y_i_t > THR_NOISE1:
                        qrs_i_raw.append(locs_temp-round(0.150*fs)+x_i_t-1)# save index of bandpass 
                        qrs_amp_raw.append(y_i_t) #save amplitude of bandpass 
                        SIG_LEV1 = 0.25*y_i_t + 0.75*SIG_LEV1 #when found with the second thres 
                    not_nois = 1
                    SIG_LEV = 0.25*pks_temp + 0.75*SIG_LEV   #when found with the second threshold              
                else:
                    not_nois = 0

        ##  find noise and QRS peaks
        if peaks[i] >= THR_SIG:
            # if a QRS candidate occurs within 360ms of the previous QRS
            # ,the algorithm determines if its T wave or QRS
            if len(qrs_c) >= 3:
                if (locs[i]-qrs_i[-1]) <= round(0.3600*fs):
                    Slope1 = np.mean(processing.calc_rr(ecg_m[locs[i]-round(0.075*fs):locs[i]])) #mean slope of the waveform at that position
                    Slope2 = np.mean(processing.calc_rr(ecg_m[qrs_i[-1]-round(0.075*fs):qrs_i[-1]])) #mean slope of previous R wave
                    if abs(Slope1) <= abs(0.5*Slope2):  # slope less then 0.5 of previous R
                        nois_c.append(peaks[i])
                        nois_i.append(locs[i])
                        skip = 1 # T wave identification
                        # adjust noise level in both filtered and
                        # MVI
                        NOISE_LEV1 = 0.125*y_i + 0.875*NOISE_LEV1
                        NOISE_LEV = 0.125*peaks[i] + 0.875*NOISE_LEV 
                    else:
                        skip = 0

            if skip == 0:  # skip is 1 when a T wave is detected       
                qrs_c.append(peaks[i])
                qrs_i.append(locs[i])
                # bandpass filter check threshold
                if y_i >= THR_SIG1:
                    if ser_back:
                        qrs_i_raw.append(x_i) # save index of bandpass 
                    else:
                        qrs_i_raw.append(locs[i]-round(0.150*fs)+x_i-1)# save index of bandpass    
                    qrs_amp_raw.append(y_i) # save amplitude of bandpass 
                    SIG_LEV1 = 0.125*y_i + 0.875*SIG_LEV1# adjust threshold for bandpass filtered sig
                # adjust Signal level
                SIG_LEV = 0.125*peaks[i] + 0.875*SIG_LEV 

        elif THR_NOISE <= peaks[i] and peaks[i] < THR_SIG:
            #adjust Noise level in filtered sig
            NOISE_LEV1 = 0.125*y_i + 0.875*NOISE_LEV1
            #adjust Noise level in MVI
            NOISE_LEV = 0.125*peaks[i] + 0.875*NOISE_LEV 

        elif peaks[i] < THR_NOISE:
            nois_c.append(peaks[i])
            nois_i.append(locs[i])
            # noise level in filtered signal
            NOISE_LEV1 = 0.125*y_i + 0.875*NOISE_LEV1 
            #adjust Noise level in MVI
            NOISE_LEV = 0.125*peaks[i] + 0.875*NOISE_LEV   

        ## adjust the threshold with SNR
        if NOISE_LEV != 0 or SIG_LEV != 0:
            THR_SIG = NOISE_LEV + 0.25*(abs(SIG_LEV - NOISE_LEV))
            THR_NOISE = 0.5*(THR_SIG)
        
        # adjust the threshold with SNR for bandpassed signal
        if NOISE_LEV1 != 0 or SIG_LEV1 != 0:
            THR_SIG1 = NOISE_LEV1 + 0.25*(abs(SIG_LEV1 - NOISE_LEV1))
            THR_NOISE1 = 0.5*(THR_SIG1)

        # take a track of thresholds of smoothed signal
        SIGL_buf.append(SIG_LEV)
        NOISL_buf.append(NOISE_LEV)
        THRS_buf.append(THR_SIG)
        # take a track of thresholds of filtered signal
        SIGL_buf1.append(SIG_LEV1)
        NOISL_buf1.append(NOISE_LEV1)
        THRS_buf1.append(THR_SIG1)

        skip = 0 #reset parameters
        not_nois = 0 #reset parameters
        ser_back = 0  #reset bandpass param   

    if gr:
        for xc in qrs_i_raw:
            ax1.axvline(x=xc, color='m', linestyle='--', linewidth=1)
            ax2.axvline(x=xc, color='m', linestyle='--', linewidth=1)
        #    ax3.axvline(x=xc, color='m', linestyle='--', linewidth=1)
            ax4.axvline(x=xc, color='m', linestyle='--', linewidth=1)
            ax5.axvline(x=xc, color='m', linestyle='--', linewidth=1)
        #ax1.plot(qrs_i_raw, qrs_amp_raw, 'o', color='y', ms=3)
        #ax2.plot(qrs_i_raw, qrs_amp_raw, 'o', color='y', ms=3)
        ax3.plot(qrs_i_raw, qrs_amp_raw, 'o', color='y', ms=3)
        #ax3.plot(locs, NOISL_buf,'--k')
        #ax3.plot(locs, SIGL_buf,'--r')
        #ax3.plot(locs,THRS_buf,'--g')
        #ax4.plot(qrs_i_raw, qrs_amp_raw, 'o', color='y', ms=3)
        #ax5.plot(qrs_i,qrs_c,'o',color='y',ms=3)
        ax6.plot(qrs_i,qrs_c,'o',color='y',ms=3)
        ax6.plot(locs, NOISL_buf,'--k',label='Noise')
        ax6.plot(locs, SIGL_buf,'--r',label='Signal Level')
        ax6.plot(locs,THRS_buf,'--g', label='Adaptive Threshold')
        plt.legend(prop={'size': 8})
        plt.tight_layout()
        plt.suptitle(subject_name)

        plt.figure()
        plt.plot(ecg,label='Raw ECG')
        plt.plot(ecg_h, label='Band-Pass Filtered')
        plt.plot(qrs_i_raw, qrs_amp_raw,'o',color='y',ms=3)
        plt.legend()

    return qrs_i_raw, qrs_amp_raw, delay, ecg_m, qrs_i, qrs_c

'''
plt.subplot(211)
plt.plot(ecg)
plt.plot(xraw,yraw, 'o',ms=3)
plt.subplot(212)
plt.plot(ecgm)
plt.plot(x,y,'o',ms=3)
plt.plot([np.mean(y)]*len(ecg))
plt.plot([.5*np.mean(y)]*len(ecg))
plt.plot([1.5*np.mean(y)]*len(ecg))
'''

def pan_plot(subject, st=0, end=3000):
    record = wfdb.rdrecord(subject)
    pan_tompkin(record.p_signal[st:end,0],record.fs,subject_name=subject)