function [] = PP_DevianceGraphs(dev, dev_smoothed, devControl, dev_min, ratio_strength, graphPairs)
%Purpose: To create a deviance graph per cell pair. Each graph will have
%         the deviance between the cell pair in black, and the control
%         deviance in red. Each graph will also denote the assembly
%         strength of the cell pair and the optimal time window for
%         prediction.

%Inputs: dev & devControl
%        dev_min (from assembly strength function)
%        ratio_strength (from assembly strength function)

%Dependencies: PP_AssemblyStrength

%Outputs:  deviance figures for specified cell pairs

%Created: 03/25/20 by Reagan Bullins
%Updated: 06/01/20 by Reagan Bullins


%% Graphing 

[~, number_of_pairs] = size(graphPairs);
pairs_to_graph = graphPairs;

min_win_pairs = zeros(1,number_of_pairs);
% Make a graph for each cell pair
for ipair = 1:number_of_pairs
   
        %Plot Dev Data, one graph for each pair
        if number_of_pairs < 5
            subplot(number_of_pairs,1,ipair)
        else
            row_num = ceil(number_of_pairs/5)
            subplot(5,row_num,ipair);
        end
        time = (0:150)
        plot(time, dev(:, pairs_to_graph(ipair)),'b')
        hold on

        plot(time, dev_smoothed(:, pairs_to_graph(ipair)),'k');
        %Find the min dev 
        %+++min_win_pairs(ipair) = find(dev(:,pairs_to_graph(ipair)) == dev_min(pairs_to_graph(ipair)));
        
            txt = (['Time Window = ' num2str(min_win_pairs(ipair)) ' ms']);
            text(1, max(devControl(:,pairs_to_graph(ipair)))+1.5, txt);
               %text(min_win_pairs(ipair) + 5, mean(dev(:,pairs_to_run(ipair))),txt)
         
            %average over all control trials in 3rd deminsion and plot
            average_control = mean(devControl,3);
            plot(time, average_control(:,pairs_to_graph(ipair)), 'r');
            
            %Label Graph Axis
            xlabel('Time (ms)')
            ylabel('Deviance')
            txt_assembly = (['Assembly Strength = ' num2str(ratio_strength(pairs_to_graph(ipair)))]);
            title({['Pair ' num2str(pairs_to_graph(ipair)) ' Deviance'], txt_assembly})
            ylim([-inf max(devControl(:,pairs_to_graph(ipair))) + 2.5])
    end
  
end