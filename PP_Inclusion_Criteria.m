% Script for cell explorer work PP_Inclusion_Criteria

%% Add paths

%Folder with spike information
    spike_info_folder = 'Recording_Specs';
%paths with result data
    path = ('C:\Users\rcbul\Documents\English Lab\');
    addpath(genpath(path));
    basepath = (['C:\Users\rcbul\Documents\English Lab\PP_RSC_Data\' spike_info_folder]);
    addpath(genpath(basepath));
    spike_info_path = [basepath 'PP_RSC_Data\' spike_info_folder];
%% Load Sorted Spikes and Raw Data
    cd(basepath);
    load('m115_191203_152410_2.spikes.cellinfo.mat'); %spike struct
    load('m115_191203_152410_2.sessionInfo.mat');
    session = sessionTemplate(pwd,'showGUI',true);
    
    %dat xml and sorted spikes
    cell_metrics = ProcessCellMetrics('session',session);
    
    cell_metrics = CellExplorer('metrics', cell_metrics);

%%
    