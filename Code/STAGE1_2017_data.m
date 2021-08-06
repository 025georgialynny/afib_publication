%% Jericho Lawson
%% Summer 2019, 2021
%% Extracting AFib Data from 2017 Dataset

% Using the mcode from the training2017 folder of the training2017.zip 
% file, the code extracts annotations from the 2017 PhysioNet/CinC 
% challenge and generates basic RR-interval information for each of 
% the subjects.

% NOTE: If you are using quantile thresholds, results of the quantile
% function will differ from results done in R/Python. This is due to
% the differences in the algorithms.

%%%% SPECIFICATIONS %%%%

WEIGHTS = [0.75, 0.25];
RUN_THRESH = [0.85, 1.15];
INTERVAL = 30;
TYPE_THRESH = "n"; % "n" for normal threshold, "q" for quantile threshold
OUT_THRESH = [0.3, 1.5];

% Directory where all files are located
SOURCE_FROM = "/Users/mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Testing/training2017/";
addpath(SOURCE_FROM); 

% Directory where new .csv files will be located
SOURCE_TO = "/Users/mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Summer 2021 Testing/Data/2017/";
SOURCE_DATA = strcat(SOURCE_TO, 'pre/');

% File of functions as object.
functs = STAGE1_2017_functions;

% Lists subjects from the 2017 PhysioNet Challenge.
fid = fopen(['RECORDS'], 'r');
RECLIST = textscan(fid, '%s');
fclose(fid);
RECORDS = RECLIST{1};

% Used to get updated list of subjects.
new_records = "";
count = 0;

% Used to gather states of each subject.
ref = readtable(strcat(SOURCE_FROM, 'REFERENCE.csv'), 'ReadVariableNames', false);

% Collects info regarding the R-R peak times and lengths; places info into
% new files.
for i = 1:length(RECORDS)
    % Collects annotations from .hea and .qrs files.
    fname = RECORDS{i};
    [tm, ecg, fs, siginfo] = rdmat(fname);
    [QRS, sign, en_thres] = qrs_detect2(ecg', 0.25, 0.6, fs); % correct
    
    % Creates starting and ending times based off RR peak information.
    times = QRS / 300;
    starts = times(1:length(QRS) - 1);
    ends = times(2:length(QRS));
    
    % Creates RR intervals.
    rr_int = diff(QRS') / 300;
    
    % Placeholders for running means, differences in RR intervals, and 
    % transition types.
    rr_mean = NaN(1, length(rr_int));
    rr_diff = NaN(1, length(rr_int));
    rr_class = strings(1, length(rr_int));
    diffs = 0;
    
    % Finds running means, differences in RR intervals, and transition types.
    for obs = 1:length(rr_int)
        % Finds if observation is outlier.
        res = functs.outlier_f(TYPE_THRESH, OUT_THRESH, rr_int(obs), rr_int);
    
        % Adds three variables to each observation if it is not an outlier.
        % Skips if it is an outlier.
        if res
            if (obs - diffs) == 1 % needed due to error with negative indices; first non-outlier
                rrm_place = 1;
            else
                rrm_place = obs - diffs - 1;
            end
            rr_mean(obs) = functs.run_mean_f(obs, diffs, WEIGHTS, rr_int(obs), rr_mean(rrm_place));
            rr_diff(obs) = functs.rr_diff_f(obs, diffs, rr_int(obs), rr_int(rrm_place));
            rr_class(obs) = functs.rr_class_f(rr_int(obs), rr_mean(obs), RUN_THRESH);
            diffs = 0;
        else
            diffs = diffs + 1;
        end
    end

    % Creates array with starting times, ending times, and RR intervals to
    % place in .csv file.
    mat = [rr_int starts.' ends.' repelem(string(ref{i, 2}), length(rr_int)).' rr_mean.' rr_diff.' rr_class.'];
    header = {'RRLength', 'Start', 'End', 'State', 'RRMean', 'RRDiff', 'RRClass'};
    if ~isempty(mat)
        mat2 = cell2table(num2cell(mat), 'VariableNames', header);
        writetable(mat2, strcat(SOURCE_DATA, string(fname), '.csv'));
        count = count + 1;
        new_records(count) = string(fname);
        
        % Data merging.
        mat2.Subject = [repelem(string(fname), length(rr_int)).'];
        if string(fname) == 'A00001'
            combined = mat2;
        else
            combined = [combined; mat2];
        end
        if string(fname) == 'A08528'
            writetable(combined, strcat(SOURCE_DATA, 'all_data_2017.csv'));
        end
    end
end

% Writes new record reference.
writematrix(new_records.', strcat(SOURCE_TO, 'ref_2017.csv'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%