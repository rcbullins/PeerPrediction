function [] = Scatter_IsoQual_by_LogLike(highest_log_value, goodUnitQual, dataset_idx1, dataset_idx2)

%Purpose: Recreate Harris' Supplement Figure 4a. Show relationship between
%         isolation distance and log likelihood of firing. Show difference
%         between two different specifications/conditions. (ie quality
%         cells versus all cells)
%Inputs:  highest_log_value (the log at the optimal time window)
%         dataset_idx1 (example: all cells)
%         dataset_idx2 (example: only quality cells)
%         unit quality
%Outputs: Scatter of relationship between isolation distance and log like.
%Dependencies: crossValidationAssemblyPrediction

%Created 7/16/20 by Reagan Bullins 
%%
plot(goodUnitQual(dataset_idx1,1), highest_log_value(dataset_idx1,1),'ob');
hold on;
plot(goodUnitQual(dataset_idx2,1), highest_log_value(dataset_idx2,1), 'or');
lc = lsline;
%lc.Color = 'r';
xlabel('Isolation Quality');
ylabel('Highest Log Likelihood');
title('Peer Predictability Dependence on Unit Quality')
legend('A','B'); % all clusters vs quality clusters

R = corrcoef(goodUnitQual(dataset_idx1,1), highest_log_value(dataset_idx1,1));
RSq = R(2)^2 
R2 = corrcoef(goodUnitQual(dataset_idx2,1), highest_log_value(dataset_idx2,1))
RSq2 = R2(2)^2

end