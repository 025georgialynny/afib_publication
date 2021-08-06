## Jericho Lawson
## Summer 2019, 2021
## Functions for Pre-Processing Program for MIT-BIH Atrial Fibrillation Project

## LIBRARIES ##
import numpy as np

## FUNCTIONS ##

# Function to determine if value is between two points.
def between_f(val, lims):
    if (val >= lims[0]) and (val <= lims[1]):
        return True
    else:
        return False

# Determines if a heartbeat is an outlier or not based off a 
# certain threshold (numerical or percentile).
def outlier_f(types, threshes, length, length_data):
    if types == "n":
        return between_f(length, threshes)
    else:
        return between_f(length, np.quantile(length_data, threshes))

def run_mean_f(ind, prev_outs, weights, curr_length, prev_mean):
    if (ind - prev_outs) == 0: # first obs
        return curr_length
    else:
        return weights[0] * prev_mean + weights[1] * curr_length

def rr_diff_f(ind, prev_outs, curr_length, prev_length):
    if (ind - prev_outs) == 0:
        return 0
    else:
        return curr_length - prev_length

def rr_class_f(val, r_mean, threshes):
    if val <= (threshes[0] * r_mean):
        return "S"
    elif val <= (threshes[1] * r_mean):
        return "Reg"
    else:
        return "L"

##############################################################################