%Assembly Script for additional plots
% Harris Supp Plots
%% Scatterplt: Iso Quality by log likelihood -HARRIS supp 4a
dataset_idx1 =  % all clusters
dataset_idx2 =  % quality clusters only
%MAKE SURE indexing is correct. goodUnitQual should be length of total
%spikes
Scatter_IsoQual_by_LogLike(highest_log_value, goodUnitQual, dataset_idx1, dataset_idx2)

%% ScatterPlot: Firing Rate by log likelihood --HARRIS Supp 4b
dataset_idx1 =  % all clusters
dataset_idx2 =  % quality clusters only

Scatter_FR_by_LogLike(spikes, optimal_win, highest_log_value, dataset_idx1, dataset_idx2)
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


