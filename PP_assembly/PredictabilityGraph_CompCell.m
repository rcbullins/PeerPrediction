function [] = PredictabilityGraph_CompCell(target_cell,log_1,log_2,optimal_win1, optimal_win2, winRange1, winRange2)
% Purpose: Make a predictability graph (log likelihood over window ranges)
%          for a singlar cell specified by the user, for two different
%          specifications (example: cell 1 while animal running vs cell 1 
%          while animal not running)
% Input:   Target Cell (Cell you want to graph)
%          log_1 (log likelihood for first specification)
%          log_2 (log likelihood for second specification)
%          optimal_win1 (optimal win for that cell in specification 1)
%          optimal_win2 (optimal win for target cell in specification 2)
%          winRange1 (winRange for 1st spec --  most likely same as winRange2)
%          winRange2
% Output:  Predictability graph with optimal time window for the target
%          cell in two different specifications/conditions.
% Dependencies: Data from assembly function
% Created: 5/16/20 by Reagan Bullins

%% Graph
% plot first data set
plot(winRange1*1000, log_1(target_cell, :), 'b');
txt2 = (['Time Window = ' num2str(optimal_win1(target_cell)*1000) ' ms']);
text(300, max(log_1(target_cell),:))-.02,txt2, 'Color','b');

% plot second data set
hold on
plot(winRange2*1000, log_2(target_cell,:), 'r');
txt2 = (['Time Window = ' num2str(optimal_win2(target_cell)*1000) ' ms']);
text(300, max(log_2(target_cell,:))-.02,txt2,'Color', 'r');

%labels
xlabel('Peer Prediction Timescale (ms)');
ylabel('Log Likelihood'); % 'Predictability (bits s-1)'
title(['Predictability for Cell:' num2str(target_cell)]); 
%xlim([0 winRange(length(winRange))*1000]);
legend('A','B');
end