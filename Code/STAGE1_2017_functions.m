%% Jericho Lawson and Georgia Smith
%% Summer 2019, 2021
%% Functions for Pre-Processing Program for MIT-BIH Atrial Fibrillation Project (2017/CinC)

classdef STAGE1_2017_functions
   methods
       
      % Determines if value is between two points.
      % I: val (point), lims (lower and upper limits)
      % O: True if point is between limits, False otherwise
      function res = between_f(obj, val, lims)
         if val >= lims(1) && val <= lims(2)
             res = true;
         else
             res = false;
         end
      end
      
      % Determines if a heartbeat is an outlier or not based off a 
      % certain threshold (numerical or percentile).
      % I: types ("n" or "q"), threshes (acceptable limits of RR intervals),
      %    length (current RR interval), length_data (all length data)
      % O: True/False depending on result of between_f call
      function res = outlier_f(obj, types, threshes, length, length_data)
          if types == "n"
              res = obj.between_f(length, threshes);
          else
              res = obj.between_f(length, quantile(length_data, threshes));
          end
      end
      
      % Finds the running mean of a certain set of values based on the current
      % value and previous running mean.
      % I: ind (current index), prev_outs (number of outliers prior to index),
      %    weights (coefficients for previous running mean and current value).
      %    curr_length (current value), prev_mean (previous running mean)
      % O: running mean for current index
      function res = run_mean_f(obj, ind, prev_outs, weights, curr_length, prev_mean)
          if (ind - prev_outs) == 1
              res = curr_length;
          else
              res = weights(1) * prev_mean + weights(2) * curr_length;
          end
      end
      
      % Finds difference in lengths between the current and previous indices.
      % I: ind (current index), prev_outs (number of outliers prior to index),
      %    curr_length (current RR interval), prev_length (previous RR interval)
      % O: 0 or difference of RR intervals
      function res = rr_diff_f(obj, ind, prev_outs, curr_length, prev_length)
          if (ind - prev_outs) == 1
              res = 0;
          else
              res = curr_length - prev_length;
          end
      end
      
      % Returns class of a certain observation based on a threshold and running mean.
      % I: val (value, RR length), r_mean (running mean), threshes (thresholds that
      %    define whether an observation is short, regular, or long)
      % O: "S", "Reg", or "L" ---> class/type of observation
      function res = rr_class_f(obj, val, r_mean, threshes)
          if val <= (threshes(1) * r_mean)
              res = "S";
          elseif val <= (threshes(2) * r_mean)
              res = "Reg";
          else
              res = "L";
          end
      end
      
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%