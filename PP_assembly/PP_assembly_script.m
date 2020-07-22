% Script for running assembly code
% Contributor: Sam
% 

%% Adding Paths

%Add Basepath for all code and data used
    basepath = ('C:\Users\rcbul\Documents\English Lab\');
%Define Recording Session Name
    %session_name = 'm115_191203_152410_n';
    session_name = 'u19_200313_155505';
%Deine DataPath that contains list of session names;
    data_path = [basepath 'PP_RSC_Data\' session_name];
% Define dataset to load/folder name
    %folder_name = '\Assembly_binLog\'; %rsc
    folder_name = '\Log_Baseline\'; %hpc
    %folder_name = '\No_IN\';
    %folder_name = '\Velocity_Baseline\';
%Define ResultPath that contains results from assembly function
    result_data_path = [data_path folder_name];
    % result_data_path = [basepath 'PP_RSC_Data\Testing\velocityAssemb']
%Add Paths
    addpath(genpath(result_data_path));
    addpath(genpath([basepath 'Code\']));
    addpath(genpath([basepath 'Sam_Code\']));
    addpath(genpath([basepath 'buzcode-dev\']));
%% Defining Specifications

%winRange = [0.002, 0.004, 0.008, 0.016, 0.032, 0.064, 0.128, 0.256, 0.512, 1.024]; %RSC assmb log
%winRange = [0.001, 0.002, 0.004, 0.008, 0.016, 0.032, 0.064, 0.128, 0.256, 0.512, 1.024]; %HPC LOG Base
 %winRange = [.001 .002:.002:.03 .032:.005:.062 .128 .256 .512 1.024] %NoIN
 
%% Run Script
% 
% for iwin = 1:length(winRange)
% [log_likelihood,weights] = CrossValidationAssemblyPrediction(spikes, 'dt', winRange(iwin), 'cells', [2]); % varargin if wanted
% save([num2str(winRange(iwin)) '_bn_VelMultiply.mat'], 'log_likelihood', 'weights','-v7.3');
% end

%% Set Graph Defaults Now
SetGraphDefaults;

%% Concat Windows & Organize Data
%NOTE: HARD CODED
[log_likelihood, weights] = Concat_Assemb_Data(winRange, result_data_path);

%% Find optimal window for every cell

[optimal_win, highest_log_value] = Find_Optimal_Window_Assemb(log_likelihood, winRange);

%% Concat EXTRA PREDICTOR
%NOTE: Hard coded CHANGE THIS CODE 
[log_velocity] = Concat_Assemb_ExtraPredic(winRange, result_data_path);

%% Find optimal window using EXTRA PREDICTOR

[optimal_win_ex, highest_log_value_ex] = Find_Optimal_Window_Assemb(log_velocity, winRange);

%% Histogram of Optimal Time Windows
bin_win_count = 1; 
winRange_graph = winRange *1000; %make in ms
bin_win_max = winRange_graph(length(winRange));
bin_win_max = bin_win_max + bin_win_count;
optimal_win_graph = optimal_win *1000; %make s to ms

[optimal_window] = PP_TimeWindow_Histogram (bin_win_count, bin_win_max, optimal_win_graph); 

%% Histogram cheat - aka bar graph with counts :)
bin_win_count = 1; 
optimal_win1 = hpc_base_optimal;
optimal_win2 = optimal_win;
winRange1 = winRangeBase;
winRange2 = winRange;

ct_per_win1 = zeros(1, length(winRange1));
ct_per_win2 = zeros(1, length(winRange2));
for iwin = 1:length(winRange1)
    ct_per_win1(1,iwin) = sum(optimal_win1(:,1) == winRange1(1,iwin));
end
for iwin = 1:length(winRange2)
    ct_per_win2(1,iwin) = sum(optimal_win2(:,1) == winRange2(1,iwin));
end

bar(1:length(winRange1), ct_per_win1, 'FaceColor', 'b', 'facealpha',.5, 'edgecolor', 'none');
hold on
bar(1:length(winRange2), ct_per_win2, 'FaceColor', 'r', 'facealpha',.5, 'edgecolor', 'none');

title('HPC vs RSC: Optimal Time Window Histogram');
ylabel('Count');
xlabel('Log Scale Time Windows (ms)');
legend('HPC', 'RSC');
set(gca,'XTickLabel', [1 2 4 8 16 32 64 128 256 512 1024])

med_opt1 = (['Median Window = ' num2str(median(optimal_win1)*1000)]);
text(.5, max(ct_per_win1), med_opt1, 'Color', 'b');
med_opt2 = (['Median Window = ' num2str(median(optimal_win2)*1000)]);
text(.5, max(ct_per_win1)-2, med_opt2, 'Color', 'r');
%% Histogram COMBO: multiple overlays of histograms for different data set -old kinda...
bin_win_count = 1; 
optimal_win1 = hpc_base_optimal;
optimal_win2 = hpc_LOG_IN_optimal;%rsc_optimal;
winRange1 = winRangeBase;
winRange2 = hpc_LOG_IN_winrange; %rsc_winRange;

winRange_graph1 = winRange1 *1000; %make in ms
bin_win_max = winRange_graph1(length(winRange1));
bin_win_max = bin_win_max + bin_win_count;
optimal_win_graph1 = optimal_win1 *1000; %make s to ms
optimal_win_graph2 = optimal_win2 *1000;


Overlay_TimeWindow_Histogram(bin_win_count, bin_win_max, optimal_win_graph1, optimal_win_graph2); 
%% predictability graph for a target cell

target_cell = 2;

plot(winRange*1000, log_likelihood(target_cell,:), 'k')
hold on
txt = (['Time Window = ' num2str(optimal_win(target_cell)*1000) ' ms']);
text(150, max(log_likelihood(target_cell,:)),txt);

xlabel('Peer Prediction Timescale (ms)');
ylabel('Log Likelihood'); % 'Predictability (bits s-1)'
title(['Predictability vs Timescale for Cell:' num2str(target_cell)]); 
xlim([0 winRange(length(winRange))*1000]);
%legend('HPC: All Cells');

%% Predictability plots for multiple dataset
index_plot = [2 4:10 12 15:40 42:45 47:60 62 63] 

log_1 = hpc_base_log;
log_2 = hpc_LOG_IN_log;
optimal_win1 = hpc_base_optimal;
optimal_win2 = hpc_LOG_IN_optimal;
winRange1 = winRangeBase;
winRange2 = hpc_LOG_IN_winrange;
for icell = 1:63
figure
target_cell = icell;

% plot first data set
plot(winRange1*1000, log_1(index_plot(target_cell), :), 'k');
txt2 = (['Time Window = ' num2str(optimal_win1(index_plot(target_cell))*1000) ' ms']);
text(300, max(log_1(index_plot(target_cell),:))-.02,txt2);

% plot second data set
hold on
plot(winRange2*1000, log_2(target_cell,:), 'r');
txt2 = (['Time Window = ' num2str(optimal_win2(target_cell)*1000) ' ms']);
text(300, max(log_2(target_cell,:))-.02,txt2,'Color', 'r');

%labels
xlabel('Peer Prediction Timescale (ms)');
ylabel('Log Likelihood'); % 'Predictability (bits s-1)'
title(['HPC: Predictability for Cell:' num2str(target_cell)]); 
%xlim([0 winRange(length(winRange))*1000]);
legend('All Cells','Only Pyramidal');
savefig([num2str(icell) '_Cell'])
end
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
    target_cell = 1;
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
plot(firing_rate(:,1), R_squared_values, 'ob');
hold on
xlabel('Firing Rate spk/s')
ylabel('R Squared')
title('All Cells: Firing Rate X R Squared')
R = corrcoef(firing_rate(:,1),R_squared_values);
R_squared_FR_Scatter = R(2)^2;
txt_fr = (['R = ' num2str(R_squared_FR_Scatter)]);
text(max(firing_rate(:,1))*.75, max(R_squared_values),txt_fr)

%% Scatter R2 X average(Abs(Weight))

plot(mean_weight_values, R_squared_values, 'ob');
hold on
xlabel('Mean Absolute Value of Peer Cell Weights')
ylabel('R Squared')
title({'All Cells: Peer Weight X R Squared'; 'HPC No Interneurons'})
R = corrcoef(mean_weight_values,R_squared_values);
R_squared_Weight_Scatter = R(2)^2;
txt_w = (['R = ' num2str(R_squared_Weight_Scatter)]);
text(max(mean_weight_values)*.75, max(R_squared_values),txt_w)
