%Script for Analyzing PP Result Data

%Inputs: Data Folder and Data Path (see first section of code to define)
%        Basepath

%Dependencies: Github
%              Peer Prediction Function Output Data

%Outputs: Assembly Strength
%         Deviance Graphs per cell pair
%         Histogram of Optimal PP Window
%         ScatterPlot of Error vs Optimal Time Window

% Created: 3/23/20 by Reagan Bullins
% Updated: 6/01/20 by Reagan Bullins          

% NOTE: Change paths to adapt to your computer
%       Go to 'EXTRA Specifications' Section to specify ... 
%           subset of data to analyze (or all)
%           pairs you want to graph (or all- but that's a lot)
%           smoothing data or not


%% Define Paths
%Folder with result data :NORMAL OR POISSON
    data_folder = 'pp_batch' ;
    %data_folder = 'pp_poisson';
%Folder with spike information
    spike_info_folder = 'spike_info';
%paths with result data
    basepath = ('C:\Users\rcbul\Documents\English Lab\');
    addpath(genpath(basepath));
    all_data_path = [basepath 'PP_RSC_Data\'];
    data_path = [basepath 'PP_RSC_Data\' data_folder];
    spike_info_path = [basepath 'PP_RSC_Data\' spike_info_folder];
%% EXTRA Specifications
  % Analyze all pairs or a subset?
        choice_analysisPairs = "All";
        %choice_analysisPairs = [2,4,5];
        %choice_analysisPairs = "WeakPairs";
        %choice_analysisPairs = "StrongPairs";
  % Graping (Deviance, Cross Corr, Rastor): Which pairs?
        choice_graphPairs = [1,51, 14];
  % Smooth Data (1) or Raw (0)?
        choice_smooth = 1;    
        %choice_smooth = 0;
  
%% Graphing Defaults

set(groot, ...
'DefaultFigureColor', 'w', ...
'DefaultAxesLineWidth', 0.5, ...
'DefaultAxesXColor', 'k', ...
'DefaultAxesYColor', 'k', ...
'DefaultAxesFontUnits', 'points', ...
'DefaultAxesFontSize', 12, ...
'DefaultAxesFontName', 'Arial', ...
'DefaultLineLineWidth', 1, ...
'DefaultTextFontUnits', 'Points', ...
'DefaultTextFontSize', 12, ...
'DefaultTextFontName', 'Arial', ...
'DefaultAxesBox', 'off', ...
'DefaultAxesTickLength', [0.02 0.025]);
 
% set the tickdirs to go out - need this specific order
set(groot, 'DefaultAxesTickDir', 'out');
set(groot, 'DefaultAxesTickDirMode', 'manual');


% To restore a property to its original MATLAB® default, use the 'remove' keyword.
% set(groot,'DefaultFigureColormap','remove')

%% Load and Concatenate Data
cd(data_path)

[dev, devControl] = PP_Load_and_Concatenate(data_path);

%% Assembly StrengthOnly 

[dev_smoothed, ratio_strength, dev_min, weak_pairs, strong_pairs] = PP_AssemblyStrength (dev, devControl, data_folder, choice_smooth);

%% Deviance Analysis 

[pairs_for_analysis, min_win_total, min_win_pairs] = PP_DevianceAnalysis (dev,dev_min choice_analysisPairs, weak_pairs, strong_pairs);

%% Deviance Graph 

PP_DevianceGraphs(dev, dev_smoothed, devControl, dev_min, ratio_strength, choice_graphPairs);

%% Histogram of Optimal Time Windows : ALL
bin_win_count = 1; %binned within 3 ms time windows -- CAN adjust
[bin_win_max,~] = size(dev); % CAN adjust 

[optimal_window] = PP_TimeWindow_Histogram (bin_win_count, bin_win_max, min_win_total); %optimal window is median

%% Histogram of Optimal Time Windows: Analysis Specifications
bin_win_count = 1; %binned within 3 ms time windows -- CAN adjust
[bin_win_max,~] = size(dev); % CAN adjust 

[optimal_window] = PP_TimeWindow_Histogram (bin_win_count, bin_win_max, min_win_pairs); %optimal window is median

%% Histogram of Assembly Strength

bin_interval = 2; 
PP_AssemblyStrength_Hist (ratio_strength, bin_interval);

%% Scatterplot of Optimal Time Windows: All

y_plot_limit = 20;
PP_OptimalWindow_ScatterPlot(min_win,dev,devControl,y_plot_limit);

%% Cross Correlograms of Pairs :Graph specifications

bin_size = .03; %also default is .01
PP_crosscorr(bin_size, choice_graphPairs,all_data_path, spike_info_path)

%% Rastor for Actual and Predictor: Graph specifications
cd(spike_info_path)
load('m115_191203_152410_2 peerPrediction_inputs.mat')
load('m115_191203_152410_2.spikes.cellinfo.mat')

cd(all_data_path)
load('pairsToRun.mat')

clear velocities
clear extraPredictors
clear position_coords

choice_sec = 105; %which second you want to plot

PP_Raster_SmoothedTrains(choice_graphPairs, choice_sec, binned_spikes, pairsToRun, spikes)

%% Place Field 



