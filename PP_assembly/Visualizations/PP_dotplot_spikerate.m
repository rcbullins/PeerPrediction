function [spikeRates] = PP_dotplot_spikerate(binned_spikes)
%Purpose: Easily visual spiking pattern across whole recording of cells.
%           Look at consistency

%Dependencies: PP_inputs.mat

%Inputs: binned_spikes (ms bins 0 or 1 : no spike or spike)

%Outputs: dot plot of spike rate (Each color a different cluster)
%         spike_rates: matrix of number of cells x minutes in recording
%                              - gives spike rate of each cell per minute

%Created: 6/10/20 by Reagan Bullins


length_rec = size(binned_spikes,3);
nbins = (1:60000:length_rec);
spike_rates = zeros(size(binned_spikes,1), length(nbins))
for cell_idx = 1:size(binned_spikes,1)
    ts = binned_spikes(cell_idx,1,:);
    ts = squeeze(ts);
    ts = ts';
    for ibin = 1:length(nbins)-1
       spike_rates(cell_idx,ibin) = sum(ts(nbins(ibin):nbins(ibin+1)));
    end
      plot(nbins, spike_rates(cell_idx,:), '.')
    hold on
end

title('Count of Spikes per minute for whole recording')
xlabel('Time (each column = 1 min)')
ylabel('Count')

end