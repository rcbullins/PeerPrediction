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

%CHANGE paths and data names if necessary for computer
%% Define Paths
%Folder with result data :NORMAL OR POISSON
    data_folder = 'pp_batch' 
    %data_folder = 'pp_poisson'
%Folder with spike information
    spike_info_folder = 'spike_info'
%paths with result data
    basepath = ('C:\English Lab\')
    addpath(genpath(basepath))
    data_path = [basepath 'PP_RSC_Data\' data_folder]
    spike_info_path = [basepath 'PP_RSC_Data\' spike_info_folder]

%% Load and Concatenate Data
cd(data_path)

[dev, devControl] = PP_Load_and_Concatenate(data_path);

%% Assembly StrengthOnly 

[ratio_strength, dev_min, dev_min_smoothed, weak_pairs, strong_pairs] = PP_AssemblyStrength (dev, devControl, data_folder);

%% Deviance Graph 
%Note: dev_min can be dev_min or dev_min_smoothed (interchangeable)
[min_win_pairs, min_win, identified_pairs] = PP_DevianceGraphs(dev, devControl, dev_min, ratio_strength, weak_pairs, strong_pairs);

%% Histogram of Optimal Time Windows : ALL
bin_win_count = 1; %binned within 3 ms time windows -- CAN adjust
[bin_win_max,~] = size(dev); % CAN adjust 

[optimal_window] = PP_TimeWindow_Histogram (bin_win_count, bin_win_max, min_win); %optimal window is median

%% Histogram of Optimal Time Windows: SPECIFIED by deviance graph function
bin_win_count = 1; %binned within 3 ms time windows -- CAN adjust
[bin_win_max,~] = size(dev); % CAN adjust 

[optimal_window] = PP_TimeWindow_Histogram (bin_win_count, bin_win_max, min_win_pairs); %optimal window is median

%% Histogram of Assembly Strength

figure
bin_win_count_ratio = 2;
nbins = (1:bin_win_count_ratio:max(ratio_strength))
histogram(ratio_strength, nbins)
[counts_per_win, edges] = histcounts(ratio_strength, nbins)
title('Assembly Strength')
xlabel('Assembly Strength')
ylabel('Frequency')
hold on
xline(median(ratio_strength), 'k','LineWidth',1.5)
txt2 = (['Median Assembly Strength = ' num2str(median(ratio_strength))]);
text(median(ratio_strength) + 5, max(counts_per_win)- 1, txt2)

%% Scatterplot of Optimal Time Windows
y_plot_limit = 20;
PP_OptimalWindow_ScatterPlot(min_win,dev,devControl,y_plot_limit);

%% Cross Correlograms of Pairs
 [~, bin_win_max] = size(position_coords);
%load spikes, and pairs to run
PP_crosscorr(identified_pairs, bin_win_max,pairsToRun, spikes)

%% Rastor for Actual and Predictor

[smoothedTrains, pair_idx] = PP_Raster_SmoothedTrains(spike_info_path, pairsToRun)

%% Place Field 



