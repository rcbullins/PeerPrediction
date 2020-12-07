function [] = getPowerSpectrum(channel, epoch)

%Purpose: To make a power spectrum graph based on given input epochs

%Dependencies: Buzcode
%              lfp mat file
%              sessionInfo

%Inputs: Channel for lfp
%        epoch start stop times [start stop; start stop; start stop]
%        session name (number)

%Outputs: Power Spectrum Graph (frequency by amplitude)

%Created: 11/18/20 by Reagan Bullins
%% Load LFP for given input channel
lfp = bz_GetLFP(channel) 

    lfp.data = lfp.data(epoch(1)*1250:epoch(2)*1250);
    lfp.timestamps = lfp.timestamps(epoch(1)*1250:epoch(2)*1250);
    lfp.interval(1) = epoch(1)
    lfp.interval(2) = epoch(2)
    
    samp_freq = lfp.samplingRate; % Set Sampling frequency
    dt = 1/samp_freq; % Sampling interval 
    rec_time = epoch(2) - epoch(1); % total recording time 
    freq_dat_1 = fft(lfp.data); % Calculate the Fourier Transform  
    pow_spec =(2*(dt^2)/rec_time)*abs(freq_dat_1); % Do the calculation for the power spectrum  
    Pow_spec2 = pow_spec(1:(length(lfp.data)/(2))+1); % Only use the beginning half of the spectrum vector  
                                  
    df = 1/max(rec_time); % Frequency resolution 
    fnq = samp_freq/2; % Nyquist frequency= half the sampling frequency.  
    freq_axis=(0:df:fnq); % Gives you the frequency axis. 
    Pow_spec2 = Pow_spec2(1:end)';
    plot(freq_axis,Pow_spec2) % Plot power spectrum 
    xlim([0 80]) % Set x-axis limits 
    xlabel('Frequency Hz')  
    ylabel ('Power') 
   
end