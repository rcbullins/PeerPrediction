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

%plot first histogram
    figure
    nbins = (1:bin_win_count:bin_win_max) 
    %histogram(min_win, nbins)
    [counts_per_win1, edges] = histcounts(optimal_win1, nbins); 
    histogram('BinEdges', edges, 'BinCounts', counts_per_win1,'FaceColor','r', 'facealpha',.5, 'edgecolor', 'none') 
    hold on
%plot second histogram
    %histogram(min_win, nbins)
    [counts_per_win2, edges] = histcounts(optimal_win2, nbins);
    histogram('BinEdges', edges, 'BinCounts', counts_per_win2, 'FaceColor', 'b', 'facealpha',.5, 'edgecolor', 'none') 
% labels
    title('Optimal Peer Prediction Time Window')
    xlabel('Optimal Time Window (ms)')
    ylabel('Count')
%identify medians for two histograms
    %xline(median(optimal_win1), '-r','LineWidth',1)
    txt1 = (['Median Window = ' num2str(median(optimal_win1)) 'ms']);
    text(20, max(counts_per_win1)- 1, txt1, 'Color', 'k')

    %xline(median(optimal_win2), '-b','LineWidth',1)
    txt2 = (['Median Window = ' num2str(median(optimal_win2)) 'ms']);
    text(20, max(counts_per_win2)- 1, txt2, 'Color', 'r')



end