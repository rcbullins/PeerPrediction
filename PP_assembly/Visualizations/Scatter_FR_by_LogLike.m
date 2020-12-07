function [] = Scatter_FR_by_LogLike(spikes, optimal_win, highest_log_value, dataset_idx1, dataset_idx2)

%Purpose: Recreate Harris' Supplement Figure 4b. Show relationship between
%         firing rate and log likelihood. Show difference
%         between two different specifications/conditions. (ie quality
%         cells versus all cells)
%Inputs:  spikes struct
%         optimal_win
%         dataset_idx1 (example: all cells)
%         dataset_idx2 (example: only quality cells)
%         highest_log_value (corresponds to optimal time window)
%Outputs: Scatter of relationship between isolation distance and log like.
%Dependencies: crossValidationAssemblyPrediction

%Created 7/16/20 by Reagan Bullins 
%%

firing_rate = zeros(length(optimal_win),1);
for icell = 1:length(optimal_win)
num_spk = length(spikes.times{icell});
length_time = spikes.times{icell}(length(spikes.times{icell})) - spikes.times{icell}(1);
firing_rate(icell,1) = num_spk/length_time;
end
plot(firing_rate(dataset_idx1,1), highest_log_value(dataset_idx1,1), 'ob');
hold on
plot(firing_rate(dataset_idx2,1), highest_log_value(dataset_idx2,1), 'or')
lc = lsline;
R = corrcoef(firing_rate(dataset_idx1,1),highest_log_value(dataset_idx1,1));
RSq = R(2)^2 
R2 = corrcoef(firing_rate(dataset_idx2,1),highest_log_value(dataset_idx2,1));
RSq2 = R2(2)^2 

xlabel('Firing Rate (spk/s)')
ylabel('Highest Log Likelihood')
title({'Peer Predictability Dependence on Firing Rate'; ['A R^2 = ' num2str(RSq) 'B R^2 = ' num2str(RSq2)]});
legend('A','B') %all vs quality or pyram vs all
end