function [optimal_win, highest_log_value] = Find_Optimal_Window_Assemb(log_likelihood, winRange)
%Purpose: Make a list of all the optimal time windows for each cell.
%Optimal time windows are defined by the window for which each cell has its
%highest log likelihood value. 
%Input: log_likelihood (num cells x winRange)
%Output: optimal_win (vector of optimal windows)
%        highest_log_value (vector of highest corresponding log values)
%Dependencies: Load_Assemb_Data function
% Created: 7/6/20 by Reagan Bullins

%%
%Initiate matrix (num cells x 1) to hold optimal time windows for each cell
num_cells = size(log_likelihood,1);
optimal_win = zeros(num_cells, 1);
%in case you want the corresponding log values that are highest
highest_log_value = zeros(num_cells, 1);

for icell = 1:num_cells
    %find highest log_likelihood out of all winRange
    highest_log_value(icell) = max(log_likelihood(icell,:));
    win_index = find(log_likelihood(icell,:) == highest_log_value(icell));
    optimal_win(icell) = winRange(win_index);
end

end