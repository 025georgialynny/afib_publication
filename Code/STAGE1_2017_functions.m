%% Jericho Lawson
%% Summer 2019, 2021
%% Functions for Pre-Processing Program for MIT-BIH Atrial Fibrillation Project

classdef STAGE1_2017_functions
   methods
      % Function to determine if value is between two points.
      function res = between_f(obj, val, lims)
         if val >= lims(1) && val <= lims(2)
             res = true;
         else
             res = false;
         end
      end
      
      % Determines if a heartbeat is an outlier or not based off a 
      % certain threshold (numerical or percentile).
      function res = outlier_f(obj, types, threshes, length, length_data)
          if types == "n"
              res = obj.between_f(length, threshes);
          else
              res = obj.between_f(length, quantile(length_data, threshes));
          end
      end
      
      function res = run_mean_f(obj, ind, prev_outs, weights, curr_length, prev_mean)
          if (ind - prev_outs) == 1
              res = curr_length;
          else
              res = weights(1) * prev_mean + weights(2) * curr_length;
          end
      end
      
      function res = rr_diff_f(obj, ind, prev_outs, curr_length, prev_length)
          if (ind - prev_outs) == 1
              res = 0;
          else
              res = curr_length - prev_length;
          end
      end
      
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