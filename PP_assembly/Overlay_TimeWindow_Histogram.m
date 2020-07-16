function [] = Overlay_TimeWindow_Histogram(bin_win_count, bin_win_max, optimal_win1, optimal_win2);

% Purpose: Plot histograms for different datasets on the same plot to make
% easy comparisons between sets of data.

%Inputs: bin_win_count (how many bins)
%        bin_win_max   (max number of bins)
%        optimal_win1  (optimal windows for dataset 1)
%        optimal_win2  (optimal windows for dataset 2)

% Dependencies:
% Output: Graph overlayed histograms showing optimal time windows for
% multiple sets of data.


hist.... ('facealpha',.5, 'edgecolor', 'none')




end