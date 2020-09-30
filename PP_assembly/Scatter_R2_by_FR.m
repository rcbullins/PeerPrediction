function [] = Scatter_R2_by_FR(optimal_win, spikes, dataset_idx1, dataset_idx2)
%Purpose: Scatter plot of R squared values (firing rate vs weight) by
%         firing rate. See if firing rate explains R squared. Also will
%         color different cell types with different colors.
%Inputs: optimal_win
%        spikes struct
%        dataset_idx1 (cell type one indexes)
%        dataset_idx2 (cell type two indexes)
%Outputs: Scatter Plot 
%Dependencies: crossValidationAssemblyPrediction
%              RSqr_FiringRate_Weights (need R squared values)

%Created 7/16/20 by Reagan Bullins 

%% 
firing_rate = zeros(length(optimal_win),1);
for icell = 1:length(optimal_win)
num_spk = length(spikes.times{icell});
length_time = spikes.times{icell}(length(spikes.times{icell})) - spikes.times{icell}(1);
firing_rate(icell,1) = num_spk/length_time;
end
interN = [1 3 11 13 14 41 46 61] %hps u19
super_inter_idx = interN %for hpc

plot(firing_rate(dataset_idx1,1), R_squared_values(dataset_idx1,1), 'ob');

hold on
plot(firing_rate(dataset_idx2,1), R_squared_values(dataset_idx2,1), 'or');

%adding pointers
% plot(firing_rate(55,1)+.6, R_squared_values(55,1), '<', 'MarkerFaceColor', 'b', 'MarkerEdgeColor','b');
% plot(firing_rate(1,1)+.6, R_squared_values(1,1), '<', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
%labels
xlabel('Target Cell Firing Rate spk/s')
ylabel('R Squared')
title({'Correlation of target cell firing rate';'to peer cell firing rate and weights'})
R = corrcoef(firing_rate(dataset_idx1,1),R_squared_values(dataset_idx1,1));
R_squared_FR_Scatter = R(2)^2;
R_p = corrcoef(firing_rate(dataset_idx2,1),R_squared_values(dataset_idx2,1));
R_squared_p = R_p(2)^2;
txt_fr = (['R = ' num2str(R_squared_FR_Scatter)]);
text(max(firing_rate(:,1))*.75, max(R_squared_values),txt_fr)
txt_p = (['R = ' num2str(R_squared_p)]);
text(max(firing_rate(:,1))*.75, max(R_squared_values)-.05,txt_p, 'Color', 'b')

legend('Dataset A','DatasetB') %pyramidal vs interneuron
end