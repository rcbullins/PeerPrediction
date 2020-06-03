function [pairs_for_analysis, min_win_total, min_win_pairs] = PP_DevianceAnalysis (dev,dev_min, analysisPairs, iweak_pairs, istrong_pairs)
%Purpose: To find the minimum optimal time window for all pairs and
%         to also find the minimum optimal time window for a specified
%         subset of pairs.

%Inputs: dev (output of peerprediction)
%        dev_min (from assembly strength function)
%        choice (specifications on what subset of data to analyze)
%        istrong_pairs & iweak_pairs

%Dependencies: PP_AssemblyStrength

%Outputs: min_win_pairs (minimum window for specified pairs) 
%         min_win_total (minimum window for all pairs)
%         pairs_for_analysis (pairs specified by choice, the index)

%Created: 6/01/20 by Reagan Bullins
%% Analysis: get min win for all windows
[~, number_of_pairs_need_minWin] = size(dev);
min_win_total = zeros(1, number_of_pairs_need_minWin);
for ipair = 1:number_of_pairs_need_minWin
    min_win_total(ipair) = find(dev(:,ipair) == dev_min(ipair));
end
    


%% Analysis: get min win for analysis specified
string_ans = isstring(analysisPairs);
if string_ans == true
    if strcmp(analysisPairs, 'WeakPairs')
        pairs_for_analysis = iweak_pairs;
    elseif strcmp(analysisPairs, 'StrongPairs')
        pairs_for_analysis = istrong_pairs;
    elseif strcmp(analysisPairs, 'All')
        pairs_for_analysis = (1:number_of_pairs_need_minWin);
    end
elseif string_ans == false
   pairs_for_analysis = analysisPairs; 
end

%find min win for specified pairs
min_win_pairs = zeros(1, length(pairs_for_analysis));
for ipair = 1:length(pairs_for_analysis)
   min_win_pairs(ipair) = find(dev(:,pairs_for_analysis(ipair)) == dev_min(pairs_for_analysis(ipair)));
end
end