% Make script with all assembly information
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
%% Make mat file with session specifications



