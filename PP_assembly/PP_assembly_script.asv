% Script for running assembly code
% 

%% Adding Paths

%Add Basepath for all code and data used
    basepath = ('C:\Users\rcbul\Documents\English Lab\');
%Define Recording Session Name
    %session_name = 'm115_191203_152410_n';
    session_name = 'u19_200313_155505'; %NUM 1 -running
    %session_name = 'u19_200310_135409'; %NUM 2
    %session_name = 'u19_200313_120452'; %NUM 3
    %session_name = 'u21_200305_153604'; %NUM 4 -runninng
    %session_name = 'u21_200309_142534'; %NUM 5
    %session_name = 'u26_200306_172032'; %NUM 6
    
%Deine DataPath that contains list of session names;
    data_path = [basepath 'PP_RSC_Data\' session_name];
% Define dataset to load/folder name
     %folder_name = '\log_fine_16_64\';
    %folder_name = '\Assembly_binLog\'; %rsc
     folder_name = '\Pulse_Epoch\';
    %folder_name = '\Log_Baseline\'; %hpc
    %folder_name = '\No_IN\';
    %folder_name = '\Velocity_Baseline\';
    %folder_name = '\DiffTimes_and_filtered\10min_nofilt_allCells\'
%Define ResultPath that contains results from assembly function
    result_data_path = [data_path folder_name];
    % result_data_path = [basepath 'PP_RSC_Data\Testing\velocityAssemb']
%Add Paths
    addpath(genpath(result_data_path));
    addpath(genpath([basepath 'Code\']));
    addpath(genpath([basepath 'Sam_Code\']));
    addpath(genpath([basepath 'buzcode-dev\']));
%% Defining Specifications
winRange = [.001 .002 .004 .008 .016:.002:.064 .128 .256 .512 1.024]; %fine_log_16_64
%winRange = [0.002, 0.004, 0.008, 0.016, 0.032, 0.064, 0.128, 0.256, 0.512, 1.024]; %RSC assmb log
%winRange = [0.001, 0.002, 0.004, 0.008, 0.016, 0.032, 0.064, 0.128, 0.256, 0.512, 1.024]; %HPC LOG Base
 %winRange = [.001 .002:.002:.03 .032:.005:.062 .128 .256 .512 1.024] %NoIN
 
%% Run Script
% 
% tic
% for iwin = 1 %:length(winRange)
% [log_likelihood,weights] = CrossValidationAssemblyPrediction(spikes, 'dt', winRange(iwin), 'epoch', [0 2400]); % varargin if wanted
% %save([num2str(winRange(iwin)) '_bn_VelMultiply.mat'], 'log_likelihood', 'weights','-v7.3');
% end
% toc

% CrossValidationAssemblyPrediction_Epochs (to give epochs)

% for ibin = 9
%     [log_likelihood, weights, log_velocity] = CrossValidationAssemblyPrediction_ExtraPredict(spikes, velocities,'dt', winRange(ibin), 'epoch', [0 2400]);
%     save([num2str(bin_list(ibin)) '_bn_asmb.mat'], 'log_likelihood', 'weights', 'log_velocity',  '-v7.3');
%     clear log_likelihood weights log_velocity
% end

%% Set Graph Defaults Now
SetGraphDefaults;

%% Concat Windows & Organize Data
%NOTE: HARD CODED
[log_likelihood, weights] = Concat_Assemb_Data(winRange, result_data_path);

% Take Away zeros -- happens because some cells do not appear til later in
% recording -CONSIDER how this could effect other graphs
%log_likelihood = log_likelihood(any(log_likelihood,2),:);

%% Find optimal window for every cell

[optimal_win, highest_log_value] = Find_Optimal_Window_Assemb(log_likelihood, winRange);

%% Concat EXTRA PREDICTOR
%NOTE: Hard coded CHANGE THIS CODE 
[log_velocity] = Concat_Assemb_ExtraPredic(winRange, result_data_path);

%% Find optimal window using EXTRA PREDICTOR

[optimal_win_ex, highest_log_value_ex] = Find_Optimal_Window_Assemb(log_velocity, winRange);

%% Histogram of Optimal Time Windows for a singular session
winRange_graph = winRange *1000; %make in ms
optimal_win_graph = optimal_win *1000; %make s to ms
%make sure axis is good for winRange-may need to alter
[optimal_window] = PP_TimeWindow_Histogram(winRange_graph, optimal_win_graph, session_name); 

%% Histogram of Optimal Time Windows for all Sessions (specified sessions)
folder_name = 'log_fine_quality_16_64';
%folder_name = 'log_fine_16_64'

path_mat_files = ['C:\Users\rcbul\Documents\English Lab\PP_RSC_Data\matFilesOverall\' folder_name];
num_mat_files = [1 2 3 5 6];

cd(path_mat_files)
Concat_Sessions_OptHist(winRange, num_mat_files)

%% Histogram Comparing two datasets (blue vs red)
optimal_win1 = 
optimal_win2 = 
winRange1 = winRange;
winRange2 = winRange; 

PP_TimeWindow_Hist_Comp(optimal_win1, optimal_win2, winRange1, winRange2)
%% Predictability graph for a target cell
target_cell = 12;

PredictabilityGraph_SingularCell(target_cell, winRange, optimal_win)
%% Predictability plots for multiple dataset
log_1 = log_likelihood;
log_2 = 
optimal_win1 = optimal_win;
optimal_win2 = ;
winRange1 = winRange;
winRange2 = winRange;
PredictabilityGraph_CompCell(target_cell,log_1,log_2,optimal_win1, optimal_win2, winRange1, winRange2)
%% Weighted Raster Plot
% DEFINE Target Cell and Second to plot on graph
    target_cell = 6;
    time_plot = 200; 
% Load spiking data
    cd(data_path)
    load([session_name '.spikes.cellinfo.mat']);
Weighted_Raster_Asmb(spikes, weights, target_cell, time_plot, optimal_win, winRange)

%% ScatterPlot : Weights by Firing Rate 
%plots the firing rate of each peer cell by the weight of that peer cell to
%the target cell
% Define a target cell
    target_cell = 55
% Load Spiking data
    cd(data_path);
    load([session_name '.spikes.cellinfo.mat']);
FiringRate_Weights_Scatter(target_cell, optimal_win, winRange, spikes, weights)

%% Distribution Graph of R values for firing rate vs weight for every cell

bin_win_count = .05;
% Load spiking data
    cd(data_path)
    load([session_name '.spikes.cellinfo.mat']);
%now gives R distr
[R_squared_values, R_values, mean_weight_values] = RSqr_FiringRate_Weights(spikes, weights, optimal_win, winRange, bin_win_count)


%% Raster sorted by weight and colored by firing rate
    figure
% DEFINE Target Cell and Second to plot on graph
    target_cell = 55;
    time_plot = 10; 
% Load spiking data
    cd(data_path)
    load([session_name '.spikes.cellinfo.mat']);

Weighted_Raster_FiringRate(spikes, weights, target_cell, time_plot, optimal_win, winRange)

%% Scatter R2 X FR
%same session (same cells) but split into different cell types
dataset_idx1 = %pyram
dataset_idx2 = %IN
Scatter_R2_by_FR(optimal_win, spikes, dataset_idx1, dataset_idx2);
%% Scatter R2 X average(Abs(Weight))
dataset_idx1 = %pyram
dataset_idx2 = %IN
Scatter_R2_by_Weight(dataset_idx1, dataset_idx2, R_squared_values, mean_weight_values)