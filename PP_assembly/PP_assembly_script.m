% Script for running assembly code
% Contributor: Sam
%
%% Adding Paths

%Folder with spike information
    spike_info_folder = 'Recording_Specs';
%paths with result data
    basepath = ('C:\Users\rcbul\Documents\English Lab\');
    addpath(genpath(basepath));
    s_path = [basepath 'Sam_Code\']; 
    spike_info_path = [basepath 'PP_RSC_Data\' spike_info_folder];
    
%% 
cd(spike_info_path);
load('m115_191203_152410_2.spikes.cellinfo.mat');
tic
[log_likelihood,weights] = CrossValidationAssemblyPrediction_Commented(spikes) % varargin if wanted
toc
cd('C:\Users\rcbul\Documents\English Lab\PP_RSC_Data\pp_assemb')
save('epoch_all.mat', 'log_likelihood', 'weights')

%% weighted raster
% Define target Cell & time window to plot
    target_cell = 1;
    time_window = 10; %what second to plot
% Get Weights of peer cells for target cell
    target_weights = (weights(1,:)); 
% Initiate indexes for length of how many peer cells there are
    idx_weights = (1:size(weights,2));
% Create matrix of indexes and weights for each peer cell
    target_idx_mat = [idx_weights;target_weights]';
% Sort the peer cells based on weights (most - to most +)
    sorted_weights = sortrows(target_idx_mat,2);
% Get index for raster plotting, flip vector so you plot most positive
% first
    raster_idx = flip(sorted_weights(:,1));

% plot target 
    h1 = subplot(2,1,1);
    title(['Target Cell: ' num2str(target_cell)])
    target_spikes = spikes.times{target_cell};
    target_x = find(target_spikes >= time_window & target_spikes < time_window+1);
        for idx_x = 1:length(target_x)
            xline(target_spikes(target_x(idx_x)));
        end
         %plot(actual_cell(actual_x),ones(length(actual_x),'.r'));
         xlim([time_window time_window+1]);
         set(gca,'XTick',[]);
         set(gca,'YTick',[]);
% make raster plot %could just plot the time window you want... this takes
% forever
    h2 = subplot(2,1,2)
    for idx_cell = 1:length(raster_idx)
       peer_cell = spikes.times{raster_idx(idx_cell)};
       predict_x = find(peer_cell >= time_window & peer_cell < time_window+1);
       y_idx = (1:length(predict_x));
       y_idx(:) = idx_cell;
       plot(peer_cell(predict_x),y_idx, '.r')
       %xline(peer_cell(predict_x))
       hold on
       xlim([time_window time_window+1]);
    end
     xlabel('Time(s)')
     ylabel('Weights')
     set(gca, 'YTick',[])
     set(h1, 'OuterPosition',[0,0.85,1,.1]);
     set(h2, 'OuterPosition',[0,.1,1,.75]);
   




