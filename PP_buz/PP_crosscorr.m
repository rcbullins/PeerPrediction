function[] = PP_crosscorr(bin_size, spikes, choice_graphPairs, pairsToRun, bin_win_max)
%Purpose: Create cross correlograms for specified pairs.

%Dependencies: Buzcode
%              Paths set in beggining of PP_Analyzing_Script
%              Spiking Data from Recording

%Inputs: choice_graphPairs (what pairs to graph)
%        bin_win_max (how long the recording is, necessary for binning)

%Output: Cross Correlogram for each pair

%Created: 06/03/20 by Reagan Bullins


%% Graph

% Get number of bins, whole length of recording split by input bin interval
%bin_num = (0:bin_interval:bin_win_max); 

for ipair = 1:length(choice_graphPairs)
    %identify which number cells are being compaired in this pair
    cell_1 = pairsToRun(choice_graphPairs(ipair), 1);
    cell_2 = pairsToRun(choice_graphPairs(ipair), 2);
    %Get the spike times for these two cells &
    %   Make an array just for them
    spiking_times{1} = spikes.times{cell_1};
    spiking_times{2} = spikes.times{cell_2};
    % Take CCG of them
    [ccg, lags] = CCG(spiking_times, [], 'binsize',bin_size)
    % CCG is 3 Dimensional, select dimension to plot
    specified_ccg = ccg(:,1,2);
    subplot(ceil(length(choice_graphPairs)/4),4,ipair);
    plot(lags, specified_ccg);
    title(['Cross Correlation: ' num2str(cell_1) ' vs ' num2str(cell_2)])
    ylabel('Count')
    xlabel('Time Lag (s)')
    xlim([-.5 .5])
end
% OLD CODE
% for ipair = 1:length(choice_graphPairs)
%     %identify which number cells are being compaired in this pair
%     cell_1 = pairsToRun(choice_graphPairs(ipair), 1);
%     cell_2 = pairsToRun(choice_graphPairs(ipair), 2);
%     %Get the spike times for these two cells
%     spikes_1 = spikes.times{cell_1};
%     spikes_2 = spikes.times{cell_2};
%     %Get counts of how many spikes in each bin
%     counts_1 = hist(spikes_1, bin_num);
%     counts_2 = hist(spikes_2, bin_num);
%     %Calculate the cross corr graph
%     [corr, lags] = xcov(counts_1, counts_2, 'coeff');
%     subplot(ceil(length(choice_graphPairs)/4),4,ipair);
%     plot(lags, corr);
%     xlim([-500 500])
%     title(['Cross Correlation: ' num2str(cell_1) ' vs ' num2str(cell_2)])
%     ylabel('Correlation')
%     xlabel('Time Lag (s)')
%     %xticklabels({'-20', '-15', '-10', '-5', '0', '5', '10', '15', '20'});
%     
% end
    
end
