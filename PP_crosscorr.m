function[] = PP_crosscorr(identified_pairs, bin_win_max, all_pairs, spikes)
% Need inputs.mat for spiking data
bin_num = (0:1:bin_win_max)

for ipair = 1:length(identified_pairs)
    cell_1 = all_pairs(identified_pairs(ipair), 1)
    cell_2 = all_pairs(identified_pairs(ipair), 2)
    
    spikes_1 = spikes.times{cell_1}
    spikes_2 = spikes.times{cell_2}
    
    counts_1 = hist(spikes_1, bin_num)
    counts_2 = hist(spikes_2, bin_num)
    
    [corr, lags] = xcov(counts_1, counts_2, 'coeff');
    subplot(ceil(length(identified_pairs)/4),4,ipair)
    plot(lags, corr);
    xlim([-10 10])
    title(['Cross Correlation: ' num2str(cell_1) ' vs ' num2str(cell_2)])
    ylabel('Correlation')
    xlabel('Time Lag (s)')
    %xticklabels({'-20', '-15', '-10', '-5', '0', '5', '10', '15', '20'});
    
end
    
end
