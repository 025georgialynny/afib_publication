%% Jericho Lawson
%% Summer 2019, 2021
%% Extracting AFib Data

% Using the mcode from the training2017 folder of the training2017.zip 
% file, the code extracts annotations from the 2017 PhysioNet/CinC 
% challenge and generates basic RR-interval information for each of 
% the subjects.

% Directory where all files are located
SOURCE_FROM = "/Users/mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Testing/training2017/";
addpath(SOURCE_FROM); 

% Directory where new .csv files will be located
SOURCE_TO = "/Users/mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Summer 2021 Testing/Data/2017/";
SOURCE_DATA = strcat(SOURCE_TO, 'pre/');

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
    RR = diff(QRS') / 300;
    
    % Creates array with starting times, ending times, and RR intervals to
    % place in .csv file.
    mat = [RR starts.' ends.' repelem(string(ref{i, 2}), length(RR)).'];
    header = {'RRLength', 'Start', 'End', 'State'};
    if ~isempty(mat)
        mat2 = cell2table(num2cell(mat), 'VariableNames', header);
        writetable(mat2, strcat(SOURCE_DATA, fname, '.csv'));
        count = count + 1;
        new_records(count) = string(RECORDS{i});
    end
end

% Writes new record reference.
writematrix(new_records.', strcat(SOURCE_TO, 'ref_2017.csv'));