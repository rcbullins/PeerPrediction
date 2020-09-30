function [] = PP_TimeWindow_Hist_Comp(optimal_win1, optimal_win2, winRange1, winRange2)
%Purpose: Creates a histogram comparing two sets of data specified by
%         optimal win 1 and 2. Each dataset will be assigned a different
%         color and specified on the legend as A or B. The medians for each
%         dataset optimal time window will also be listed.
%Inputs: optimal_win1
%        optimal_win2 
%        winRange1
%        winRange2
%Outputs: Histogram of optimal peer prediciton time windows comparing two
%         different datasets
%Dependencies: find optimal window of assemb function
%              output of CrossValidationAssemblyPrediction
%Created 7/16/20 by Reagan Bullins

%% Find number of cells that have each time windows as its optimal
ct_per_win1 = zeros(1, length(winRange1));
ct_per_win2 = zeros(1, length(winRange2));
for iwin = 1:length(winRange1)
    ct_per_win1(1,iwin) = sum(optimal_win1(:,1) == winRange1(1,iwin));
end
for iwin = 1:length(winRange2)
    ct_per_win2(1,iwin) = sum(optimal_win2(:,1) == winRange2(1,iwin));
end
%% Graph
bar(1:length(winRange1), ct_per_win1, 'FaceColor', 'b', 'facealpha',.5, 'edgecolor', 'none');
hold on
bar(1:length(winRange2), ct_per_win2, 'FaceColor', 'r', 'facealpha',.5, 'edgecolor', 'none');

title('Optimal Time Window Histogram');
ylabel('Count');
xlabel('Log Scale Time Windows (ms)');
legend('A', 'B');
set(gca,'XTickLabel', [1 2 4 8 16 32 64 128 256 512 1024])

med_opt1 = (['Median Window = ' num2str(median(optimal_win1)*1000)]);
text(.5, max(ct_per_win1), med_opt1, 'Color', 'b');
med_opt2 = (['Median Window = ' num2str(median(optimal_win2)*1000)]);
text(.5, max(ct_per_win1)-2, med_opt2, 'Color', 'r');

end