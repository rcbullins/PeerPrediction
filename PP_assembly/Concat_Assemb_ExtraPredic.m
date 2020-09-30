function [log_extraPredictor] = Concat_Assemb_ExtraPredic(winRange, result_data_path)
% Purpose: Load and concatenate data from assembly function 
% Input: winRange, result_data_path
% Output: log_liklihood extra predictor = (num cell x winRange)
% Dependencies: Data from assembly function
% Created: 7/10/20 by Reagan Bullins

%% 
%go to data folder
cd(result_data_path);
% Concat data to make a log likelihood matrix (cell x winRange)
%first get num cells
    load([num2str(winRange(1)) '_bn_Vel2.mat'])
    num_cells = size(log_likelihood,1);
% initiate matrix and array
    log_likelihood_all = zeros(num_cells,length(winRange));
% For every window
for iwin = 1:length(winRange)
    % load the data for that given window
    load([num2str(winRange(iwin)) '_bn_Vel2.mat'])
    % add data to concat matrix and array
    log_likelihood_all(:,iwin) = log_velocity;
end
    log_extraPredictor = log_likelihood_all;
   
    
end