function [spk_whole_epoch, spk_center_epoch] = spks_per_epoch(spikes, pulseEpochs)
% Purpose: Count how many times each cell fires within an epoch (pulses
%          given experimentally or something similar). Count how many times the cell
%          fires the whole epoch and also just the center of the epoch.
% Inputs:  spikes (spikes struct, need the times)
%          pulseEpochs(start and stop times of epochs)
% Dependencies: currently epochs are known to be 300 ms, so center of epoch
%               is epoch start time + 100 ms and epoch stop time - 100 ms.
%               --can change this if needed.
% Outputs: spk_whole_epoch (gives matrix number of cells x number of
%                           epochs) gives count for number of spikes
%          spk_center_epoch (same as whole, except only number of spikes in
%                            center of epoch, ie 100 ms after and before
%                            epoch begin and end times)
% Created: 7/21/20 by Reagan Bullins


%%
% Create matrix num cells x num epochs, for the whole epoch
spk_whole_epoch = zeros(length(spikes.times), length(pulseEpochs));
% Create another matrix num cells x num epochs, for the center of the epoch
spk_center_epoch = zeros(length(spikes.times), length(pulseEpochs));

% for every cell ...
for icell = 1:length(spikes.times)
    %check every epoch to see if icell spikes & how many times
    for iepoch = 1:length(pulseEpochs)
        %for whole epoch, how many times does this cell spike
        spk_indexes = find(spikes.times{icell}(:) >= pulseEpochs(iepoch,1) ...
                        & spikes.times{icell}(:) <= pulseEpochs(iepoch,2));
        spk_whole_epoch(icell, iepoch) = length(spk_indexes);
        %for center of epoch, how many times does this cell spike
        spk_center_indexes = find(spikes.times{icell}(:) >= pulseEpochs(iepoch,1)+.1 ...
                               & spikes.times{icell}(:) <= pulseEpochs(iepoch,2)-.1);
        spk_center_epoch(icell, iepoch) = length(spk_center_indexes);
    end
end
%% figure out consistency of firing

        
        whole_epoch_spk_percentage = zeros (length(spikes.times),1);
        center_epoch_spk_percentage = zeros(length(spikes.times),1);
        %make spike either 1 or 0 for yes there are spikes or no there isnt
        spk_whole_iteration = spk_whole_epoch;
        spk_center_iteration = spk_center_epoch;
        spk_whole_iteration(spk_whole_iteration >= 1) = 1;
        spk_center_iteration(spk_center_iteration >= 1) = 1;
        
        whole_epoch_spk_percentage = zeros(length(spikes.times),1);
        center_epoch_spk_percentage = zeros(length(spikes.times),1);
        
        whole_epoch_spk_percentage = sum(spk_whole_iteration,2)/size(pulseEpochs,1)
        center_epoch_spk_percentage = sum(spk_center_iteration,2)/size(pulseEpochs,1)
  
        %visualize percentage of epochs in
        nbins = (0:.1:1)
        [N, edges] = histcounts(whole_epoch_spk_percentage, nbins)
        histogram('BinEdges', edges, 'BinCounts', N);
        title('Whole Epoch Participation by cell')
        xlabel('Percentage of Epoch Participation')
        ylabel('Cell Count')
        
        %visualize percentage of epochs in
        nbins = (0:.1:1)
        [N, edges] = histcounts(center_epoch_spk_percentage, nbins)
        histogram('BinEdges', edges, 'BinCounts', N);
        title('Center Epoch Participation by cell')
        xlabel('Percentage of Epoch Participation')
        ylabel('Cell Count')
%% heat map
imagesc(spk_whole_epoch)

imagesc(spk_center_epoch)

%% wanna see by epoch number
ct_per_whole_epoch = sum(spk_whole_iteration,1);
ct_per_center_epoch = sum(spk_center_iteration,1);

bar(1:843, ct_per_whole_epoch(1,:));
title('Cell participation by epoch')
ylabel('Number of Cells')
xlabel('Epoch number');
%% only interneurons
          interN = [1 3 11 13 14 41 46 61]
        for iIN = 1:length(interN)
            spikesIN{iIN} = spikes.times{interN(iIN)};
        end
        whole_epoch_spk_percentage = zeros (length(spikesIN),1);
        center_epoch_spk_percentage = zeros(length(spikesIN),1);
        %make spike either 1 or 0 for yes there are spikes or no there isnt
        spk_whole_iteration = spk_whole_epoch(interN,:);
        spk_center_iteration = spk_center_epoch(interN,:);
        spk_whole_iteration(spk_whole_iteration >= 1) = 1;
        spk_center_iteration(spk_center_iteration >= 1) = 1;
        
        whole_epoch_spk_percentage = zeros(length(spikesIN),1);
        center_epoch_spk_percentage = zeros(length(spikesIN),1);
        
        whole_epoch_spk_percentage = sum(spk_whole_iteration,2)/size(pulseEpochs,1)
        center_epoch_spk_percentage = sum(spk_center_iteration,2)/size(pulseEpochs,1)
  
        %visualize percentage of epochs in
        nbins = (0:.1:1)
        [N, edges] = histcounts(whole_epoch_spk_percentage, nbins)
        histogram('BinEdges', edges, 'BinCounts', N);
        title({'Whole Epoch Participation by cell';'Interneurons'})
        xlabel('Percentage of Epoch Participation')
        ylabel('Cell Count')
        
        %visualize percentage of epochs in
        nbins = (0:.1:1)
        [N, edges] = histcounts(center_epoch_spk_percentage, nbins)
        histogram('BinEdges', edges, 'BinCounts', N);
        title({'Center Epoch Participation by cell';'Interneurons'})
        xlabel('Percentage of Epoch Participation')
        ylabel('Cell Count')
        

%% find how many times each cell fires all together in all epochs
       total_epoch_length = length(pulseEpochs) * .3; %in seconds 
       sum_center = sum(spk_center_epoch, 2);
       sum_whole = sum(spk_whole_epoch, 2);
       
       fr_center = sum_center/total_epoch_length; %firing rate spk/s
       fr_whole = sum_whole/total_epoch_length;
       
       %QUESTION: only include firing rate above a certain number???
%% Visualize Spike Counts
 % for whole epoch
 bar(1:length(spikes.times), sum_center(:,1));
 title('Cell spike count sum for center of epochs')
 xlabel('Cell Number')
 ylabel('Number of Spikes')
 
 figure
 bar(1:length(spikes.times), sum_whole(:,1));
 title('Cell spike count sum for whole epoch')
 xlabel('Cell Number')
 ylabel('Number of Spikes')
 
 %% see what's up in the baseline during this time
 
 base_whole_epoch = zeros(length(spikes.times), length(pulseEpochs));
 base_center_epoch = zeros(length(spikes.times), length(pulseEpochs));
 
 for icell = 1:length(spikes.times)
    %check every epoch to see if icell spikes & how many times
    for iepoch = 1:length(pulseEpochs)
        %for whole epoch, how many times does this cell spike
        base_indexes = find(spikes.times{icell}(:) >= pulseEpochs(iepoch,1)- 2400 ...
                        & spikes.times{icell}(:) <= pulseEpochs(iepoch,2)- 2400);
        base_whole_epoch(icell, iepoch) = length(base_indexes);
        %for center of epoch, how many times does this cell spike
        base_center_indexes = find(spikes.times{icell}(:) >= pulseEpochs(iepoch,1)+.1 -2400 ...
                               & spikes.times{icell}(:) <= pulseEpochs(iepoch,2)-.1 -2400);
        base_center_epoch(icell, iepoch) = length(base_center_indexes);
    end
 end
 
%% Visualize baseline (-40 min from epochs)
       total_epoch_length = length(pulseEpochs) * .3; %in seconds 
       base_sum_center = sum(base_center_epoch, 2);
       base_sum_whole = sum(base_whole_epoch, 2);
       
       base_fr_center = base_sum_center/total_epoch_length; %firing rate spk/s
       base_fr_whole = base_sum_whole/total_epoch_length;
       

bar(1:length(spikes.times), base_sum_center(:,1));
 title('Baseline: Cell spike count for center of epochs')
 xlabel('Cell Number')
 ylabel('Number of Spikes')
 
 figure
 bar(1:length(spikes.times), base_sum_whole(:,1));
 title('Baseline: Cell spike count for whole epochs')
 xlabel('Cell Number')
 ylabel('Number of Spikes')
 

end


