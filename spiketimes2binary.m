function [binned_spikes] = spiketimes2binary(params, spikeTimes)
  % Creat 1 ms bins over length of recording
    bin_num = (0:.001:params.length_recording);
%     maxT = max(cellfun(@max,spikes.times))
%     bin_num = (0:.001:maxT);
%     
    %Make matrix of zeros(bins x spike times)
    binned_spikes = zeros(length(spikeTimes),2, length(bin_num)-1); %subtract 1 from length from binning
    %For each cell, bin spikes by 1 ms, turn into binary
    for cell_num = 1: length(spikeTimes)
        [counts,bin_edges] = histcounts(spikeTimes{cell_num}, bin_num);
        %If the count in a bin is greater than 0, the cell fired, make it 1
        counts(counts > 0) = 1;
        binned_spikes(cell_num,1,:) = counts;   
        binned_spikes(cell_num,2,:) = counts;
    end
end 