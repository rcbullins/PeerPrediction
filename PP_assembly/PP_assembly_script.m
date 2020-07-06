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


%% weighted raster
% Define target Cell & second of recording to plot
    target_cell = 1;
    time_plot = 10; %what second to plot
% Get Weights of peer cells for target cell
    target_weights = (weights{target_cell}(target_cell+1,:)); 
    idx_weights = (1:num_cells);
% Take out nan value (given target cell idx)
    idx_nan = find(isnan(target_weights));
    target_weights(idx_nan) = [];
    idx_weights(idx_nan) = [];
    zero_index = find(target_weights < 0, 1, 'first');
% Create matrix of indexes and corresponding weights for each peer cell
    target_idx_mat = [idx_weights;target_weights]';
% Sort the peer cells based on weights (most - to most +)
    sorted_weights = sortrows(target_idx_mat,2);
% Add in row of zeros where weights change negative to positive: graphing
% purposes
    zero_index = find(sorted_weights(:,2) > 0, 1, 'first');
    sorted_weights_wZ = zeros(num_cells,2);
    sorted_weights_wZ(:,1) = [sorted_weights(1:zero_index-1,1)' 0 sorted_weights(zero_index:length(sorted_weights),1)']';
    sorted_weights_wZ(:,2) = [sorted_weights(1:zero_index-1,2)' 0 sorted_weights(zero_index:length(sorted_weights),2)']';
% Get index for raster plotting
    raster_idx = sorted_weights_wZ(:,1);

%% plot target cell
    h1 = subplot(2,1,1);
    title(['Target Cell: ' num2str(target_cell)])
    target_spikes = spikes.times{target_cell};
    target_x = find(target_spikes >= time_plot & target_spikes < time_plot+1);
        for idx_x = 1:length(target_x)
            xline(target_spikes(target_x(idx_x)));
        end
         xlim([time_plot time_plot+1]);
         set(gca,'XTick',[]);
         set(gca,'YTick',[]);
%% make raster plot
    h2 = subplot(2,1,2)
    % Assigning color to row (yes important, yes kinda lengthy)
    %RGB = rgb('dark red','magenta', 'rose', 'purplish', 'burple', 'true blue')
   
    % for each cell, plot a raster
    for idx_cell = 1:length(raster_idx)
       % if raster_idx is ZERO, this is a marker on the graph
       if raster_idx(idx_cell) == 0
           yline(idx_cell, '--');
       else
       %current cell spikes (positive to negative weight)
       peer_cell = spikes.times{raster_idx(idx_cell)};
       % only plot given window to plot
       predict_x = find(peer_cell >= time_plot & peer_cell < time_plot+1);
       % get y-axis established
       y_idx = (1:length(predict_x));
       % which y value to plot this cell on
       y_idx(:) = idx_cell;
       plot(peer_cell(predict_x),y_idx, '.r')
       %xline(peer_cell(predict_x))
       hold on
       xlim([time_plot time_plot+1]);
    end
    end 
     xlabel('Time(s)')
     ylabel('Weights')
     set(gca, 'YTick',[])
     set(h1, 'OuterPosition',[0,0.85,1,.1]);
     set(h2, 'OuterPosition',[0,.1,1,.75]);
    



