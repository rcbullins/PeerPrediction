function [] = FiringRate_Weights_Scatter(target_cell, optimal_win, winRange, spikes, weights)

%Purpose: Creates a scatter plot comparing firing rate(x axis) by weight of
%         peer cell (y axis). Additionally, a fitted line of the scatter
%         plot is depicted to show the trend of the data.
%Inputs: target_cell (defined by the user)
%        optimal_win (optimal window given by highest log likelihood for a cell)
%        winRange (window bin sizes)
%        spikes struct (really only need spikes.times)
%        weights (output of crossvalidationassemblyPrediction)
%Outputs: Scatterplot of weights by firing rate
%Dependencies: find optimal window of assemb function
%              output of CrossValidationAssemblyPrediction
%Created 7/16/20 by Reagan Bullins

%%
% Find optimal window for target_cell
    tar_opt_win = optimal_win(target_cell,1);
% Find which index in winRange this window is
    win_idx = find(winRange == tar_opt_win)
% Get Weights of peer cells for target cell
    target_weights = (weights{win_idx}(target_cell+1,:)); 
    num_cells = size(target_weights, 2)
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
% Get index for raster plotting
    raster_idx = sorted_weights(:,1);
% Get average firing rate for each of these cells
    firing_rate = zeros(length(raster_idx), 1);
    max_time = max(cellfun(@max,spikes.times));
    for icell = 1:length(raster_idx)
        num_spikes = length(spikes.times{raster_idx(icell)});
        firing_rate(icell, 1) = num_spikes/max_time;
    end
%Plot firing rate by weight
f = fit(firing_rate, sorted_weights(:,2), 'poly1')
plot(f,firing_rate, sorted_weights(:,2), 'o');
xlabel('Firing Rate (spikes/s)');
ylabel('Weights of Peer Cells');
title(['Target Cell: ' num2str(target_cell)]);
xlim([0 max(firing_rate(:,1))]);
legend('Peer Cell', 'Fitted');
end