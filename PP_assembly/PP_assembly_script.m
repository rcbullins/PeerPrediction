% Script for running assembly code
% Contributor: Sam
% 

%% Adding Paths

%Add Basepath for all code and data used
    basepath = ('C:\Users\rcbul\Documents\English Lab\');
%Define Recording Session Name
    session_name = 'm115_191203_152410_n';
    %session_name = 
%Deine DataPath that contains list of session names;
    data_path = [basepath 'PP_RSC_Data\' session_name];
%Define ResultPath that contains results from assembly function
    result_data_path = [data_path '\Assembly_binLog\'];
    % result_data_path = [basepath 'PP_RSC_Data\Testing\velocityAssemb']
%Add Paths
    addpath(genpath(result_data_path));
    addpath(genpath([basepath 'Code\']));
    addpath(genpath([basepath 'Sam_Code\']));
    addpath(genpath([basepath 'buzcode-dev\']));
%% Defining Specifications

winRange = [0.002, 0.004, 0.008, 0.016, 0.032, 0.064, 0.128, 0.256, 0.512, 1.024]; %RSC assmb log
%winRange = [0.001, 0.002, 0.004, 0.008, 0.016, 0.032, 0.064, 0.128, 0.256, 0.512, 1.024]; %HPC LOG Base
 %winRange = [.001 .002:.002:.2]
 
%% Run Script
% cd(data_path);
% load([session_name '.spikes.cellinfo.mat']);
% tic
% for iwin = 1:length(winRange)
% [log_likelihood,weights, log_velocity] = CrossValidationAssemblyPrediction_ExtraPredict(spikes, velocities, 'dt', winRange(iwin), 'cells', [2]); % varargin if wanted
%  save([num2str(winRange(iwin)) '_bn_VelMultiply.mat'], 'log_likelihood', 'weights', 'log_velocity','-v7.3');
% end
% toc

for iwin = 1:length(winRange)
[log_likelihood,weights] = CrossValidationAssemblyPrediction(spikes, 'dt', winRange(iwin), 'cells', [2]); % varargin if wanted
save([num2str(winRange(iwin)) '_bn_VelMultiply.mat'], 'log_likelihood', 'weights','-v7.3');
end
% redefining spikes to use
% spikes_IN = {};
% pyram_cells = [2 4:10 12 15:40 42:45 47:60 62:63]
% 
% for i = 1:length(pyram_cells)
% spikes_IN{i} = spikes.times{pyram_cells(i)};
% end
%% Set Graph Defaults Now
SetGraphDefaults;

%% Concat Windows & Organize Data
%NOTE: HARD CODED
[log_likelihood, weights] = Concat_Assemb_Data(winRange, result_data_path);

%% Concat Extra Predictors
%NOTE: Hard coded CHANGE THIS CODE 
[log_velocity] = Concat_Assemb_ExtraPredic(winRange, result_data_path);
%% Find optimal window for every cell

[optimal_win, highest_log_value] = Find_Optimal_Window_Assemb(log_likelihood, winRange);

%% Find optimal window using EXTRA PREDICTOR

[optimal_win_ex, highest_log_value_ex] = Find_Optimal_Window_Assemb(log_velocity, winRange);

%% Histogram of Optimal Time Windows
bin_win_count = 1; 
winRange_graph = winRange *1000; %make in ms
bin_win_max = winRange_graph(length(winRange));
bin_win_max = bin_win_max + bin_win_count;
optimal_win_graph = optimal_win *1000; %make s to ms

[optimal_window] = PP_TimeWindow_Histogram (bin_win_count, bin_win_max, optimal_win_graph); 

%% Histogram COMBO: multiple overlays of histograms for different data set
bin_win_count = 1; 
winRange_graph = winRange *1000; %make in ms
bin_win_max = winRange_graph(length(winRange));
bin_win_max = bin_win_max + bin_win_count;
optimal_win_graph1 = optimal_win1 *1000; %make s to ms
optimal_win_graph2 = optimal_win2 *1000;

[optimal_window] = Overlay_TimeWindow_Histogram(bin_win_count, bin_win_max, optimal_win_graph1, optimal_win_graph2); 
%% predictability graph for a target cell
target_cell = 2;

plot(winRange*1000, log_likelihood(target_cell,:), 'k')
hold on
txt = (['Time Window = ' num2str(optimal_win(target_cell)*1000) ' ms']);
text(300, max(log_likelihood(target_cell,:)),txt);

% plot(winRange*1000, log_velocity(target_cell,:), 'r');
% txt2 = (['Time Window = ' num2str(optimal_win_ex(target_cell)*1000) ' ms']);
% text(300, max(log_likelihood(target_cell,:))-.02,txt2);

xlabel('Peer Prediction Timescale (ms)');
ylabel('Log Likelihood'); % 'Predictability (bits s-1)'
title(['Predictability vs Timescale for Cell:' num2str(target_cell)]); 
%xlim([0 winRange(length(winRange))*1000]);
legend('Peer Activity Alone', 'Peer Activity + Velocity');


%% Weighted Raster Plot
% DEFINE Target Cell and Second to plot on graph
    target_cell = 1;
    time_plot = 200; 
% Load spiking data
    cd(data_path)
    load([session_name '.spikes.cellinfo.mat']);
Weighted_Raster_Asmb(spikes, weights, target_cell, time_plot, optimal_win, winRange)

%% ScatterPlot : Weights by Firing Rate
% Define a target cell
    target_cell = 10;
% Load Spiking data
    cd(data_path);
    load([session_name '.spikes.cellinfo.mat']);
FiringRate_Weights_Scatter(target_cell, optimal_win, winRange, spikes, weights)

%% Distribution Graph of R^squared values for firing rate vs weight for every cell

bin_win_count = .05;
% Load spiking data
    cd(data_path)
    load([session_name '.spikes.cellinfo.mat']);

[R_squared_values, mean_weight_values] = RSqr_FiringRate_Weights(spikes, weights, optimal_win, winRange, bin_win_count)

%% Raster sorted by weight and colored by firing rate
    figure
% DEFINE Target Cell and Second to plot on graph
    target_cell = 4;
    time_plot = 5000; 
% Load spiking data
    cd(data_path)
    load([session_name '.spikes.cellinfo.mat']);

Weighted_Raster_FiringRate(spikes, weights, target_cell, time_plot, optimal_win, winRange)

%% Scatter R2 X FR
firing_rate = zeros(length(optimal_win),1);
for icell = 1:length(optimal_win)
num_spk = length(spikes.times{icell});
length_time = spikes.times{icell}(length(spikes.times{icell})) - spikes.times{icell}(1);
firing_rate(icell,1) = num_spk/length_time;
end

%% Scatter R2 X average(Abs(Weight))

plot(mean_weight_values, R_squared_values, 'ob');
xlabel('Mean Absolute Value of Peer Cell Weights')
ylabel('R Squared')
title('All Cells: Peer Weight X R Squared')
R = corrcoef(mean_weight_values,R_squared_values);
R_squared_Weight_Scatter = R(2)^2;
