% Script for running assembly code
% Contributor: Sam
% 

%% Adding Paths

%Add Basepath for all code and data used
    basepath = ('C:\Users\rcbul\Documents\English Lab\');
    addpath(genpath(basepath));
%Define Recording Session Name
    session_name = 'm115_191203_152410_n';
%Deine DataPath that contains list of session names;
    data_path = [basepath 'PP_RSC_Data\' session_name];
%Define ResultPath that contains results from assembly function
    result_data_path = [data_path '\Assembly\'];

%% Defining Specifications

winRange = [0.002, 0.004, 0.008, 0.016, 0.032, 0.064, 0.128, 0.256, 0.512, 1.024];
    
%% Run Script
% cd(data_path);
% load([session_name '.spikes.cellinfo.mat']);
% tic
% [log_likelihood,weights] = CrossValidationAssemblyPrediction_Commented(spikes) % varargin if wanted
% toc
% cd('C:\Users\rcbul\Documents\English Lab\PP_RSC_Data\pp_assemb')
% save('epoch_all.mat', 'log_likelihood', 'weights')

%% Set Graph Defaults Now
SetGraphDefaults;

%% Concat Windows & Organize Data

[log_likelihood, weights] = Concat_Assemb_Data(winRange, result_data_path);

%% Find optimal window for every cell

[optimal_win, highest_log_value] = Find_Optimal_Window_Assemb(log_likelihood, winRange);

%% Histogram of Optimal Time Windows
bin_win_count = 1; 
winRange_graph = winRange *1000; %make in ms
bin_win_max = winRange_graph(length(winRange));
bin_win_max = bin_win_max + bin_win_count;
optimal_win_graph = optimal_win *1000; %make s to ms

[optimal_window] = PP_TimeWindow_Histogram (bin_win_count, bin_win_max, optimal_win_graph); 

%% predictability graph for a target cell
target_cell = 2;

plot(winRange*1000, log_likelihood(target_cell,:))
hold on
txt = (['Time Window = ' num2str(optimal_win(target_cell)*1000) ' ms']);
text(300, max(log_likelihood(target_cell,:)),txt);
xlabel('Peer Prediction Timescale (ms)');
ylabel('Log Likelihood'); % 'Predictability (bits s-1)'
title(['Predictability vs Timescale for Cell:' num2str(target_cell)]); 
xlim([0 winRange(length(winRange))*1000]);


%% Weighted Raster Plot
% DEFINE Target Cell and Second to plot on graph
    target_cell = 20;
    time_plot = 200; 
% Load spiking data
    cd(data_path)
    load([session_name '.spikes.cellinfo.mat']);
Weighted_Raster_Asmb(spikes, weights, target_cell, time_plot, optimal_win, winRange)
