function [] = getPowerSpectrum_Advanced(channel, epochRun, epochNoRun)
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

lfp = bz_GetLFP(channel) % hard coded 57 for 155505
for i = 1:length(epochRun)
lfpRun.data = lfp.data(epochRun(i,1)*1250:epochRun(i,2)*1250);
lfpRun.timestamps = lfp.timestamps(epochRun(i,1)*1250:epochRun(i,2)*1250);
lfpRun.interval(i,1) = epochRun(i,1)
lfpRun.interval(i,2) = epochRun(i,2)
end
lfp = bz_GetLFP(channel) % hard coded 57 for 155505
for i = 1:length(epochNoRun)
lfpNoRun.data = lfp.data(epochNoRun(i,1)*1250:epochNoRun(i,2)*1250);
lfpNoRun.timestamps = lfp.timestamps(epochNoRun(i,1)*1250:epochNoRun(i,2)*1250);
lfpNoRun.interval(i,1) = epochNoRun(i,1)
lfpNoRun.interval(i,2) = epochNoRun(i,2)
end
%% 
samp_freq = lfp.samplingRate; % Set Sampling frequency
dt = 1/samp_freq; % Sampling interval 

% First for running
rec_time_run = zeros(1, length(epochRun))
for itime = 1:length(epochRun)
rec_time_run(1,itime) = epochRun(itime,2) - epochRun(itime,1); % total recording time 
end
rec_time_run_tot = sum(rec_time_run);

freq_run = fft(lfpRun.data); % Calculate the Fourier Transform  
pow_spec_run =(2*(dt^2)/rec_time_run_tot)*abs(freq_run); % Do the calculation for the power spectrum  
Pow_spec2_run = pow_spec_run(1:(length(lfpRun.data)/(2))+1); % Only use the beginning half of the spectrum vector  
  
% For not running
rec_time_no = zeros(1, length(epochNoRun))
for itime = 1:length(epochNoRun)
rec_time_no(1,itime) = epochNoRun(itime,2) - epochNoRun(itime,1); % total recording time 
end
rec_time_no_tot = sum(rec_time_no);

freq_noRun = fft(lfpNoRun.data); % Calculate the Fourier Transform  
pow_spec_no =(2*(dt^2)/rec_time_no_tot)*abs(freq_noRun); % Do the calculation for the power spectrum  
Pow_spec2_no = pow_spec_no(1:(length(lfpNoRun.data)/(2))+1); % Only use the beginning half of the spectrum vector  
  
% Plotting
df = 1/max(rec_time_run_tot); % Frequency resolution 
fnq = samp_freq/2; % Nyquist frequency= half the sampling frequency.  
freq_axis=(0:df:fnq); % Gives you the frequency axis. 
Pow_spec2_run = Pow_spec2_run(1:end-1)';
Pow_spec2_no = Pow_spec2_no(1:end-1)';
% WORK HERE - add in freq axis as x-axis
 % Plot power spectrum 
lh = plot(Pow_spec2_no)
hold on
lh.Color = [1,0,0,0.1]
hold on
lf = plot(Pow_spec2_run)
lf.Color = [0,0,1,0.1]
%xlim([0 80]) % Set x-axis limits 
xlabel('Frequency Hz')  
ylabel ('Power') 
legend('No Movement','Movement')