function [] = Scatter_R2_by_Weight(dataset_idx1, dataset_idx2, R_squared_values, mean_weight_values)
%Purpose: Scatter plot of R squared values (firing rate vs weight) by
%         weight. See if weight explains R squared. Also will
%         color different cell types with different colors.
%Inputs: R_squared values
%        mean_weight_values
%        dataset_idx1 (cell type one indexes)
%        dataset_idx2 (cell type two indexes)
%Outputs: Scatter Plot 
%Dependencies: crossValidationAssemblyPrediction
%              RSqr_FiringRate_Weights (need R squared values)

%Created 7/16/20 by Reagan Bullins 

%% 
plot(mean_weight_values(dataset_idx1,1), R_squared_values(dataset_idx1,1), 'ob');
hold on
plot(mean_weight_values(dataset_idx2,1), R_squared_values(dataset_idx2,1), 'or');
%specify pointers
% plot(mean_weight_values(55,1)+.05, R_squared_values(55,1), '<', 'MarkerFaceColor', 'b', 'MarkerEdgeColor','b');
% plot(mean_weight_values(1,1)+.05, R_squared_values(1,1), '<', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');

xlabel('Mean Absolute Value of Peer Cell Weights')
ylabel('R Squared')
title({'Correlation of target cell weights';'to peer cell firing rate and weights'})

R = corrcoef(mean_weight_values(dataset_idx1,1),R_squared_values(dataset_idx1,1));
R_squared_Weight_Scatter = R(2)^2;
R_p = corrcoef(mean_weight_values(dataset_idx2,1),R_squared_values(dataset_idx2,1));
R_squared_p = R_p(2)^2;
txt_w = (['R^2 = ' num2str(R_squared_Weight_Scatter)]);
text(max(mean_weight_values)*.75, max(R_squared_values),txt_w)
txt_p = (['R^2 = ' num2str(R_squared_p)]);
text(max(mean_weight_values)*.75, max(R_squared_values)-.05,txt_p, 'Color','b')
legend('A','B') % Pyram vs IN
end