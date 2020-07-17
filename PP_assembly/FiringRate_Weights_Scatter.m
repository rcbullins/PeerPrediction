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
% Get Weights of peer cells for target cell (+1 is because first row is
% constant)
    target_weights = (weights{win_idx}(target_cell+1,:)); 
    num_cells = size(target_weights, 2)
    idx_weights = (1:num_cells);
% Take out nan value (given target cell idx)
    idx_nan = find(isnan(target_weights));
    target_weights(idx_nan) = [];
    idx_weights(idx_nan) = [];
% Create matrix of indexes and corresponding weights for each peer cell
    target_idx_mat = [idx_weights;abs(target_weights)]';
% Sort the peer cells based on weights (most - to most +)
    sorted_weights = sortrows(target_idx_mat,2);
% Get index of cells
    raster_idx = sorted_weights(:,1);
% Get average firing rate for each of these cells
    firing_rate = zeros(length(raster_idx), 1);
    for icell = 1:length(raster_idx)
        %get length of time from first spike to last
        firing_time = spikes.times{raster_idx(icell)}(length(spikes.times{raster_idx(icell)})) - spikes.times{raster_idx(icell)}(1);
        %get number of spikes for this cell
        num_spikes = length(spikes.times{raster_idx(icell)});
        %get average spiking rate for this cell
        firing_rate(icell, 1) = num_spikes/firing_time;
    end
    firing_time_target = spikes.times{target_cell}(length(spikes.times{target_cell})) - spikes.times{target_cell}(1);
    num_spikes_target = length(spikes.times{target_cell});
    FR_tar = num_spikes_target/firing_time_target;
%Get R Squared Value and fitted line
R = corrcoef(firing_rate,sorted_weights(:,2));
R_squared = R(2)^2 
% 
% [fitted, S, mu] = polyfit(firing_rate, sorted_weights(:,2), 1);
% f = polyval(fitted, firing_rate, S, mu);
% plot(firing_rate, sorted_weights(:,2), 'o', firing_rate, f, '-r');

plot(firing_rate, sorted_weights(:,2), 'ob');
hold on
lc = lsline;
lc.Color = 'r';
yline(.25, '--');
xlabel('Firing Rate (spikes/s)');
ylabel('Abs Weights of Peer Cells');
title({['Target Cell: ' num2str(target_cell)]; ['R^2 = ' num2str(R_squared) ' & FR = ' num2str(FR_tar) ' spk/s']});
xlim([0 max(firing_rate(:,1))]);
legend('Peer Cell', 'LSR', '.25 weight');
end