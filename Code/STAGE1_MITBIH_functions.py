## Jericho Lawson and Georgia Smith
## Summer 2019, 2021
## Functions for Pre-Processing Program for MIT-BIH Atrial Fibrillation Project 
## (MIT-BIH)

#### LIBRARIES ####
import numpy as np

#### FUNCTIONS ####

# Determines if value is between two points.
# I: val (point), lims (lower and upper limits)
# O: True if point is between limits, False otherwise
def between_f(val, lims):
    if (val >= lims[0]) and (val <= lims[1]):
        return True
    else:
        return False

# Determines if a heartbeat is an outlier or not based off a 
# certain threshold (numerical or percentile).
# I: types ("n" or "q"), threshes (acceptable limits of RR intervals),
#    length (current RR interval), length_data (all length data)
# O: True/False depending on result of between_f call
def outlier_f(types, threshes, length, length_data):
    if types == "n":
        return between_f(length, threshes)
    else:
        return between_f(length, np.quantile(length_data, threshes))

# Finds the running mean of a certain set of values based on the current
# value and previous running mean.
# I: ind (current index), prev_outs (number of outliers prior to index),
#    weights (coefficients for previous running mean and current value).
#    curr_length (current value), prev_mean (previous running mean)
# O: running mean for current index
def run_mean_f(ind, prev_outs, weights, curr_length, prev_mean):
    if (ind - prev_outs) == 0: # first obs
        return curr_length
    else:
        return weights[0] * prev_mean + weights[1] * curr_length

# Finds difference in lengths between the current and previous indices.
# I: ind (current index), prev_outs (number of outliers prior to index),
#    curr_length (current RR interval), prev_length (previous RR interval)
# O: 0 or difference of RR intervals
def rr_diff_f(ind, prev_outs, curr_length, prev_length):
    if (ind - prev_outs) == 0:
        return 0
    else:
        return curr_length - prev_length

# Returns class of a certain observation based on a threshold and running mean.
# I: val (value, RR length), r_mean (running mean), threshes (thresholds that
#    define whether an observation is short, regular, or long)
# O: "S", "Reg", or "L" ---> class/type of observation
def rr_class_f(val, r_mean, threshes):
    if val <= (threshes[0] * r_mean):
        return "S"
    elif val <= (threshes[1] * r_mean):
        return "Reg"
    else:
        return "L"

##############################################################################