function [] = PP_AssemblyStrength_Hist (ratio_strength,bin_interval)

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

%Outputs: Histogram of Assembly Strength

%Created: 06/02/20 by Reagan Bullins


figure
bin_win_count_ratio = bin_interval;
nbins = (1:bin_win_count_ratio:max(ratio_strength))
[counts_per_win, edges] = histcounts(ratio_strength, nbins)
histogram('BinEdges', edges, 'BinCounts', counts_per_win);

title('Assembly Strength')
xlabel('Assembly Strength')
ylabel('Frequency')
hold on
xline(median(ratio_strength), 'k','LineWidth',1.5)
txt2 = (['Median Assembly Strength = ' num2str(median(ratio_strength))]);
text(median(ratio_strength) + 5, max(counts_per_win)- 1, txt2)
end