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
% Edited: 6/01/20 by Reagan Bullins          

% NOTE: Change paths to adapt to your computer
%       Go to 'EXTRA Specifications' Section to specify inclusion/exclusion
%       criteria.

%% Define Paths
%Folder with result data :NORMAL OR POISSON
    data_folder = 'pp_batch' 
    %data_folder = 'pp_poisson'
%Folder with spike information
    spike_info_folder = 'spike_info'
%paths with result data
    basepath = ('C:\Users\rcbul\Documents\English Lab\PP_RSC_Data\');
    addpath(genpath(basepath))
    data_path = [basepath data_folder]
    spike_info_path = [basepath spike_info_folder]
%% EXTRA Specifications
  % Analyze all pairs or a subset?
        %choice_analysisPairs = "All"
        %choice_analysisPairs = [2,4,5];
        choice_analysisPairs = "WeakPairs"
        %choice_analysisPairs = "StrongPairs"
  % Graping (Deviance and Cross Corr): Which pairs?
        choice_graphPairs = [2,4,5];
        
%% Load and Concatenate Data
cd(data_path)

[dev, devControl] = PP_Load_and_Concatenate(data_path);

%% Assembly StrengthOnly 

[ratio_strength, dev_min, dev_min_smoothed, weak_pairs, strong_pairs] = PP_AssemblyStrength (dev, devControl, data_folder);

%% Deviance Analysis 
[pairs_for_analysis, min_win_total, min_win_pairs] = PP_DevianceAnalysis (dev,dev_min, choice_analysisPairs, weak_pairs, strong_pairs);

%% Deviance Graph 
%Note: dev_min can be dev_min or dev_min_smoothed (interchangeable)
PP_DevianceGraphs(dev, devControl, dev_min, ratio_strength, choice_graphPairs);

%% Histogram of Optimal Time Windows : ALL
bin_win_count = 1; %binned within 3 ms time windows -- CAN adjust
[bin_win_max,~] = size(dev); % CAN adjust 

[optimal_window] = PP_TimeWindow_Histogram (bin_win_count, bin_win_max, min_win_total); %optimal window is median

%% Histogram of Optimal Time Windows: Analysis Specifications
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

%% Scatterplot of Optimal Time Windows: All
y_plot_limit = 20;
PP_OptimalWindow_ScatterPlot(min_win,dev,devControl,y_plot_limit);

%% Cross Correlograms of Pairs :Graph specifications
 [~, bin_win_max] = size(position_coords);
 cd(basepath)
 load(pairsToRun.mat)
%load spikes, and pairs to run
PP_crosscorr(choice.graphPairs, bin_win_max ,pairsToRun, spikes)

%% Rastor for Actual and Predictor: Graph specifications
cd(basepath)
 load(pairsToRun.mat)
[smoothedTrains, pair_idx] = PP_Raster_SmoothedTrains(spike_info_path, pairsToRun)

%% Place Field 



