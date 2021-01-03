% Script to create a mat file to run data over ARC
%% Adding Paths
%Add Basepath for folder with code and data folders
    basepath = (''); 
%Define Recording Session Name
    session_name = ''
%Deine DataPath that contains list of session names;
    data_folder = ''
    code_folder = ''
    data_path = [basepath '\' data_folder '\' session_name];
%Add Paths
    addpath(genpath([basepath code_folder '\']));
    addpath(genpath([basepath 'buzcode-dev\']));

%% Make mat file with Globus Information
cd(data_path)
% Load spike struct
load([session_name '.spikes.cellinfo.mat'])
% Define time windows to run modeling code over, we suggest the following
bin_list = [.001 .002 .004 .008 .016:.002:.064 .128 .256 .512 1.024];
% Define a folder you want to save your data to on globus
results_folder = ''
% Define what epoch of time you want to run the modeling code over
start_time = ; % in seconds
stop_time = ; % in seconds
epochRun = [start_time stop_time]

save([session_name '_globus.mat', 'spikes','session_name','results_folder','epochRun','bin_list']);
