function [optimal_window] = PP_TimeWindow_Histogram (winRange, min_win)
%Purpose: To make a histogram of the optimal time windows for all cell
%         pairs OR more specifically by the ones denoted from the previous
%         functions -- PP_DevianceGraphs. The histogram will have the mode
%         of time windows denoted with a black line and the median denoted
%         by a red line.

%Dependencies: PP_DevianceGraphs
%              PP_AssemblyStrength

%Inputs: bin_win_count (how many windows per bin)
%        bin_win_max (maximum time window to graph)
%        min_win (the minimum dev corresponding to optimal time window)

%Outputs: optiml_window (the optimal time window for predicition)
%            --CURRENTLY: the bin with max count. Median? Mode? 

%Created: 3/31/20 by Reagan Bullins

%%

ct_per_win = zeros(1, length(winRange));

for iwin = 1:length(winRange)
    ct_per_win(1,iwin) = sum(min_win(:,1) == winRange(1,iwin));
end
bar(1:length(winRange), ct_per_win, 'FaceColor', 'b', 'facealpha',.5, 'edgecolor', 'none');
hold on
title('Optimal Time Window Histogram');
ylabel('Count');
xlabel('Log Scale Time Windows (ms)');
set(gca,'XTickLabel', [1 2 4 8 16 32 64 128 256 512 1024])
med_opt1 = (['Median Window = ' num2str(median(min_win))]);
text(.5, max(ct_per_win), med_opt1, 'Color', 'r');

optimal_window = median(min_win);
%% old iteration
% prompt = 'Histogram graph on log scale? Yes OR No: '
% choice_log = input(prompt,'s');
% 
% figure
% nbins = (1:bin_win_count:bin_win_max) 
% %histogram(min_win, nbins)
% [counts_per_win, edges] = histcounts(min_win, nbins)
% histogram('BinEdges', edges, 'BinCounts', counts_per_win) %normalize??
% title('Optimal Peer Prediction Time Window')
% xlabel('Optimal Time Window (ms)')
% ylabel('Count')
% hold on
% xline(median(min_win), 'k','LineWidth',1.5)
% txt2 = (['Median Window = ' num2str(median(min_win)) 'ms']);
% text(20, max(counts_per_win)- 1, txt2)
% 
% optimal_window = median(min_win);
% 
% if strcmp(choice_log, 'Yes')
%     set(gca, 'Yscale', 'log')
%     ylabel('Log10(Frequency)')
% end
end