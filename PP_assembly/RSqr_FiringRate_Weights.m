function [R_squared, R_values, avg_abs_weight] = RSqr_FiringRate_Weights(spikes, weights, optimal_win, winRange, bin_win_count)

%Purpose: Creates a histogram of R squared values for every cell against
%its peer cells, for firing rate vs weight. 
%Inputs: spikes (spikes struct)
%        weights (given by CrossValidationAssemblyPrediction function)
%        optimal_win (optimal window per cell)
%Outputs: Histogram of RSquared Values of abs weights vs firing rate
%Dependencies: crossValidationAssemblyPrediction

%Created 7/17/20 by Reagan Bullins 
%%

R_squared = zeros(length(optimal_win),1);
R_values = zeros(length(optimal_win),1);
avg_abs_weight = zeros(length(optimal_win),1);

for target_cell = 1:length(optimal_win)
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
    
%Get R Squared Value and fitted line
R = corrcoef(firing_rate,sorted_weights(:,2));
R_values(target_cell,1) = R(2);
R_squared(target_cell,1) = R(2)^2
avg_abs_weight(target_cell,1) = mean(abs(target_weights));
end

%%make histogram
figure
nbins = (-.5:bin_win_count:max(R_values(:,1)));
%histogram(min_win, nbins)
[counts_per_win, edges] = histcounts(R_values(:,1), nbins)
histogram('BinEdges', edges, 'BinCounts', counts_per_win, 'FaceColor', 'b', 'facealpha',.5)
title('Histogram of R: |Weights| by Firing Rate')
xlabel('R')
ylabel('Count')
hold on
xline(median(R_values(:,1)), 'k','LineWidth',1.5)
txt2 = (['Median R: ' num2str(median(R_values))])
text(.4, max(counts_per_win)- 1, txt2)




end

