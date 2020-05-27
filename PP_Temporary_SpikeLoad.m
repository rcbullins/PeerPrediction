%% temporary loading in spikes for analysis
    params.nChans               = sessionInfo.nChannels;
    params.sampFreq             = sessionInfo.rates.wideband;
    params.juxtaspikes_noTM     = 1;
    params.length_recording = max(spikes.times{:})
    %% Load in Spikes and lfp
    spikes = bz_GetSpikes;
    %spikes = bz_GetSpikes('sortingMethod','kilosort');
    channel1_lfp = bz_GetLFP(0); % for time of recording reference
    params.length_recording = max(channel1_lfp.timestamps);
    %%
 [binned_spikes] = spiketimes2binary(params, spikes.times)