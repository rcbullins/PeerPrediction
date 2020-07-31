%Assembly Script for additional plots
% Harris Supp Plots
%% Scatterplt: Iso Quality by log likelihood -HARRIS supp 4a

plot(goodUnitQual(:,1), highest_log_value(:,1),'ob');
hold on;
plot(goodUnitQual(supergoodUnitQualIdx(:,1)), highest_log_value(supergoodUnitQualIdx(:,1)), 'or');
lc = lsline;
%lc.Color = 'r';
xlabel('Isolation Quality');
ylabel('Highest Log Likelihood');
title('HPC: Peer Predictability Dependence on Unit Quality')
legend('All Clusters','Quality Clusters');

R = corrcoef(goodUnitQual(supergoodUnitQualIdx(:,1)), highest_log_value(supergoodUnitQualIdx(:,1)));
R_squared = R(2)^2 
R2 = corrcoef(goodUnitQual(:,1), highest_log_value(:,1))
RSq2 = R2(2)^2
text(1,1,'R^2 = 0.3968', 'Color', 'r')
text(1,1,'R^2 = 0.4372', 'Color', 'k')
%%%%%%%%%%%%%%%%%%Same but for only pyram

plot(goodUnitQual(pyram_idx,1), highest_log_value(pyram_idx,1),'ob');
hold on;
super_idx_pyram = [2 4 5 6 10:20 22 23 25:32 34 35];
plot(goodUnitQual(supergoodUnitQualIdx(super_idx_pyram,1)), highest_log_value(supergoodUnitQualIdx(super_idx_pyram,1)), 'or');
lc2 = lsline;
%lc.Color = 'r';
xlabel('Isolation Quality');
ylabel('Highest Log Likelihood');
title({'HPC: Peer Predictability Dependence on Unit Quality';'Target Cells: Only Pyramidal'})
legend('All Clusters','Quality Clusters');

R = corrcoef(goodUnitQual(pyram_idx,1), highest_log_value(pyram_idx,1));
R_squared = R(2)^2 
R2 = corrcoef(goodUnitQual(supergoodUnitQualIdx(super_idx_pyram,1)), highest_log_value(supergoodUnitQualIdx(super_idx_pyram,1)))
RSq2 = R2(2)^2
text(1,1,'R^2 = 0.0441', 'Color', 'b')
text(1,1,'R^2 = 0.0183', 'Color', 'r')

%%%% extra
qual = goodUnitQual(supergoodUnitQualIdx(super_idx_pyram,1))
high_log = highest_log_value(supergoodUnitQualIdx(super_idx_pyram,1))
tbl = table(qual, high_log);
[z,p,k] = filtm(tbl, 'qual~high_log')

B1 = [ones(size(lc1.XData(:))), lc1.XData(:)]\lc1.YData(:);
slope_all = B1(2)

B2 = [ones(size(lc2.XData(:))), lc2.XData(:)]\lc2.YData(:);
slope_sup = B2(2)
anova(lm2, 'summary')
%% ScatterPlot: Firing Rate by log likelihood --HARRIS Supp 4b

firing_rate = zeros(length(optimal_win),1);
for icell = 1:length(optimal_win)
num_spk = length(spikes.times{icell});
length_time = spikes.times{icell}(length(spikes.times{icell})) - spikes.times{icell}(1);
firing_rate(icell,1) = num_spk/length_time;
end
plot(firing_rate(:,1), highest_log_value(:,1), 'ob');
hold on
%plot(firing_rate(supergoodUnitQualIdx,1), highest_log_value(supergoodUnitQualIdx,1), 'or')
plot(firing_rate(pyram_super_idx,1), highest_log_value(pyram_super_idx,1), 'or')
lc = lsline;
%R = corrcoef(firing_rate(supergoodUnitQualIdx,1),highest_log_value(supergoodUnitQualIdx,1));
R = corrcoef(firing_rate(pyram_super_idx,1),highest_log_value(pyram_super_idx,1));
R_squared = R(2)^2 
%lc.Color = 'r';
xlabel('Firing Rate (spk/s)')
ylabel('Highest Log Likelihood')
title({'HPC: Peer Predictability Dependence on Firing Rate'; ['R^2 = ' num2str(R_squared)]});
legend('Pryamidal Clusters >= 20 Isolation Quality')

%% Scatter of Log Likelihood for each time of recording 
idx_time_pyram = [2 4 5 6 10:20 22 23 25:32 34 35]

min10_x = zeros(length(min10_highest_log(idx_time_pyram,1)),1);
min10_x(:,1) = 1;
min20_x = zeros(length(min20_highest_log(idx_time_pyram,1)),1);
min20_x(:,1) = 2;
min30_x = zeros(length(min30_highest_log(idx_time_pyram,1)),1)
min30_x(:,1) = 3;
min40_x = zeros(length(min40_highest_log(idx_time_pyram,1)),1)
min40_x(:,1) = 4;

plot(min10_x, min10_highest_log(idx_time_pyram,1), '.r')
hold on
plot(min20_x, min20_highest_log(idx_time_pyram,1), '.r')
plot(min30_x, min30_highest_log(idx_time_pyram,1), '.r')
plot(min40_x, min40_highest_log(idx_time_pyram,1), '.r')
xticks([1 2 3 4])
xticklabels({'10','20','30','40'})
xlim([0 5])
xlabel('Recording Length (min)')
ylabel('Log Likelihood')
title('HPC: Peer Predictability Dependence on Recording Length')


%% Scatter for recording length vs optimal window

min10_x = zeros(length(min10_optimal_win),1);
min10_x(:,1) = 1;
min20_x = zeros(length(min20_optimal_win),1);
min20_x(:,1) = 2;
min30_x = zeros(length(min30_optimal_win),1)
min30_x(:,1) = 3;
min40_x = zeros(length(min40_optimal_win),1)
min40_x(:,1) = 4;

plot(min10_x, min10_optimal_win(:,1), '.r')
hold on
plot(min20_x, min20_optimal_win(:,1), '.r')
plot(min30_x, min30_optimal_win(:,1), '.r')
plot(min40_x, min40_optimal_win(:,1), '.r')
xticks([1 2 3 4])
xticklabels({'10','20','30','40'})
xlim([0 5])
xlabel('Recording Length (min)')
ylabel('Optimal Window')
title('HPC: Optimal Window Dependence on Recording Length')

%% Bar graph for different recording lengths
idx_time_pyram = [2 4 5 6 10:20 22 23 25:32 34 35]
ct_per_win1 = zeros(1, length(winRange));
ct_per_win2 = zeros(1, length(winRange));
ct_per_win3 = zeros(1, length(winRange));
ct_per_win4 = zeros(1, length(winRange));
for iwin = 1:length(winRange)
    ct_per_win1(1,iwin) = sum(min10_optimal_win(idx_time_pyram,1) == winRange(1,iwin));
    ct_per_win2(1,iwin) = sum(min20_optimal_win(idx_time_pyram,1) == winRange(1,iwin));
    ct_per_win3(1,iwin) = sum(min30_optimal_win(idx_time_pyram,1) == winRange(1,iwin));
    ct_per_win4(1,iwin) = sum(min40_optimal_win(idx_time_pyram,1) == winRange(1,iwin));
end

bars_time = zeros(length(winRange),4);
bars_time(:,1) = ct_per_win1;
bars_time(:,2) = ct_per_win2;
bars_time(:,3) = ct_per_win3;
bars_time(:,4) = ct_per_win4;

b = bar(bars_time)
legend('10 min','20 min', '30 min', '40 min')
xticklabels({'1','2','4','8','16','32','64','128','256', '512','1024'})
xlabel('Optimal Time Window (ms)')
ylabel('Count')
title({'Optimal Time Window dependence on Recording Length'; 'Good Isolation Quality: Only Pyramidal'})
hold on
b(1).FaceColor = [0 .447 .741];
b(2).FaceColor = [.85 .325 .098];
b(3).FaceColor = [.929 .694 .125];
b(4).FaceColor = [.494 .184 .556];

txt1 = (['Median Window = ' num2str(median(min10_optimal_win(:,1)))]);
txt2 = (['Median Window = ' num2str(median(min20_optimal_win(:,1)))]);
txt3 = (['Median Window = ' num2str(median(min30_optimal_win(:,1)))]);
txt4 = (['Median Window = ' num2str(median(min40_optimal_win(:,1)))]);
text(.5, max(ct_per_win1)-1.5, txt1, 'Color', [0 .447 .741]);
text(.5, max(ct_per_win1)-3, txt2, 'Color', [.85 .325 .098]);
text(.5, max(ct_per_win1)-4.5, txt3, 'Color', [.929 .694 .125]);
text(.5, max(ct_per_win1)-6, txt4, 'Color', [.494 .184 .556]);

%% Histogram for different windows
bin_win_count = 1; 
optimal_win1 = min40_optimal_win;
optimal_win2 = min30_optimal_win;
winRange1 = winRange;
winRange2 = winRange;

ct_per_win1 = zeros(1, length(winRange1));
ct_per_win2 = zeros(1, length(winRange2));
for iwin = 1:length(winRange1)
    ct_per_win1(1,iwin) = sum(optimal_win1(:,1) == winRange1(1,iwin));
end
for iwin = 1:length(winRange2)
    ct_per_win2(1,iwin) = sum(optimal_win2(:,1) == winRange2(1,iwin));
end

bar(1:length(winRange1), ct_per_win1, 'FaceColor', 'b', 'facealpha',.5, 'edgecolor', 'none');
hold on
bar(1:length(winRange2), ct_per_win2, 'FaceColor', 'r', 'facealpha',.5, 'edgecolor', 'none');

title({'HPC: Optimal Time Window Histogram';,'No Inclusion Factor'});
ylabel('Count');
xlabel('Log Scale Time Windows (ms)');
legend('40 min', '30 min');
set(gca,'XTickLabel', [1 2 4 8 16 32 64 128 256 512 1024])

med_opt1 = (['Median Window = ' num2str(median(optimal_win1)*1000)]);
text(.5, max(ct_per_win1), med_opt1, 'Color', 'b');
med_opt2 = (['Median Window = ' num2str(median(optimal_win2)*1000)]);
text(.5, max(ct_per_win1)-2, med_opt2, 'Color', 'r');
