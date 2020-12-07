function [log_likelihood, weights] = Concat_Assemb_Data(winRange, result_data_path)
% Purpose: Load and concatenate data from assembly function 
% Input: winRange, result_data_path
% Output: log_liklihood = (num cell x winRange)
%         weights = {winRange}(num cell+1 x num cell)
% Dependencies: Data from assembly function
% Created: 7/6/20 by Reagan Bullins

%% 
%go to data folder
cd(result_data_path);
% Concat data to make a log likelihood matrix (cell x winRange)
%first get num cells
    load([num2str(winRange(1)) '_bn_asmb.mat'])
    num_cells = size(log_likelihood,1);
% initiate matrix and array
    log_likelihood_all = zeros(num_cells,length(winRange));
    weights_all = {};
% For every window
for iwin = 1:length(winRange)
    % load the data for that given window
    load([num2str(winRange(iwin)) '_bn_asmb.mat'])
    % add data to concat matrix and array
    log_likelihood_all(:,iwin) = log_likelihood;
    weights_all{iwin} = weights;
end
    log_likelihood = log_likelihood_all;
    weights = weights_all;
    
end