function [] = PredictabilityGraph_SingularCell(target_cell, winRange, optimal_win, log_likelihood)
% Purpose: Make a predictability graph (log likelihood over window ranges)
%          for a singlar cell specified by the user. Lists the optimal
%          window for this cell on the graph.
% Input:   Target Cell (Cell you want to graph)
%          winRange (window Ranges used)
%          optimal_win (optimal win for that cell)
% Output:  Predictability graph with optimal time window for the target
%          cell.
% Dependencies: Data from assembly function
% Created: 5/16/20 by Reagan Bullins
%% Graph Log Likelihood
plot(winRange*1000, log_likelihood(target_cell,:), 'r')

hold on
txt = (['Time Window = ' num2str(optimal_win(target_cell)*1000) ' ms']);
text(150, max(log_likelihood(target_cell,:)),txt);

xlabel('Peer Prediction Timescale (ms)');
ylabel('Log Likelihood'); % 'Predictability (bits s-1)'
title(['Predictability vs Timescale for Cell:' num2str(target_cell)]); 
xlim([0 winRange(length(winRange))*1000]);
end