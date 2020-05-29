function [min_win_pairs, min_win_total, pairs_to_run] = PP_DevianceGraphs(dev, devControl, dev_min, ratio_strength, istrong_pairs, iweak_pairs)
%Purpose: To create a deviance graph per cell pair. Each graph will have
%         the deviance between the cell pair in black, and the control
%         deviance in red. Each graph will also denote the assembly
%         strength of the cell pair and the optimal time window for
%         prediction.
s
%Inputs: dev & devControl
%        dev_min (from assembly strength function)
%        ratio_strength (from assembly strength function)

%Dependencies: PP_AssemblyStrength

%Outputs: min_win (optimal pp window per cell pair)
%         deviance figures for specified cell pairs

%Created: 3/25/20 by Reagan Bullins
   

prompt = 'Do you want to graph All cell pairs or Specify? All OR Specify: '
choice_all = input(prompt, 's');


if strcmp(choice_all, 'All')
    [~,number_of_pairs] = size(dev);
elseif strcmp(choice_all, 'Specify')
    prompt_strength = 'Do you want to compare by assembly strength? Yes OR No: '
    choice_strength = input(prompt_strength, 's');
    if strcmp(choice_strength, 'Yes')
        prompt_str_wk = 'Which assembly strength do you want to run? Strong OR Weak: '
        choice_str_wk = input(prompt_str_wk, 's');
        if strcmp(choice_str_wk, 'Strong')
            pairs_to_run = istrong_pairs;
        elseif strcmp(choice_str_wk, 'Weak')
            pairs_to_run = iweak_pairs;
        end
    elseif strcmp(choice_strength, 'No') 
        prompt2 = 'Which cell pairs do you want to run? State Answer as such [#,#,#,..] : '
        pairs_to_run = input(prompt2)
    end
    number_of_pairs = length(pairs_to_run);
end

%% get min win for all windows
[~, number_of_pairs_need_minWin] = size(dev);
min_win_total = zeros(1, number_of_pairs_need_minWin);
for ipair = 1:number_of_pairs_need_minWin
    min_win_total(ipair) = find(dev(:,ipair) == dev_min(ipair));
end
    

%%
min_win_pairs = zeros(1,number_of_pairs);
% Make a graph for each cell pair

for ipair = 1:number_of_pairs
    if strcmp(choice_all, 'All')
         %Plot Dev Data, one graph for each pair
            row_num = ceil(number_of_pairs/5);
            subplot(5,row_num,ipair);
            time = (0:150)
         plot(time, dev(:,ipair), 'b');
         hold on
         %Find the min dev 
            min_win_pairs(ipair) = find(dev(:,ipair) == dev_min(ipair));
            txt = (['Time Window = ' num2str(min_win_pairs(ipair)) ' ms']);
            text(1, max(devControl(:,ipair))+1.5, txt);
             %text(min_win_pairs(ipair) + 5, dev(min_win_pairs(ipair)),txt)
           
            %average over all control trials in 3rd deminsion and plot
            average_control = mean(devControl,3);
            plot(time, average_control(:,ipair), 'r');
            %Label Graph Axis
            xlabel('Time (ms)')
            ylabel('Deviance')
            txt_assembly = (['Assembly Strength = ' num2str(ratio_strength(ipair))]);
            title({['Pair ' num2str(ipair) ' Deviance'], txt_assembly})
            ylim([-inf max(devControl(:,ipair))+ 2.5])
            
    elseif strcmp(choice_all, 'Specify')
        %Plot Dev Data, one graph for each pair
        if number_of_pairs < 5
            subplot(number_of_pairs,1,ipair)
        else
            row_num = ceil(number_of_pairs/5)
            subplot(5,row_num,ipair);
        end
        time = (0:150)
        plot(time, dev(:, pairs_to_run(ipair)),'b')
        hold on
        %Find the min dev 
        min_win_pairs(ipair) = find(dev(:,pairs_to_run(ipair)) == dev_min(pairs_to_run(ipair)));
        
            txt = (['Time Window = ' num2str(min_win_pairs(ipair)) ' ms']);
            text(1, max(devControl(:,pairs_to_run(ipair)))+1.5, txt);
               %text(min_win_pairs(ipair) + 5, mean(dev(:,pairs_to_run(ipair))),txt)
         
            %average over all control trials in 3rd deminsion and plot
            average_control = mean(devControl,3);
            plot(time, average_control(:,pairs_to_run(ipair)), 'r')
            %Label Graph Axis
            xlabel('Time (ms)')
            ylabel('Deviance')
           txt_assembly = (['Assembly Strength = ' num2str(ratio_strength(pairs_to_run(ipair)))]);
            title({['Pair ' num2str(pairs_to_run(ipair)) ' Deviance'], txt_assembly})
            ylim([-inf max(devControl(:,pairs_to_run(ipair))) + 2.5])
    end
  
end
end 