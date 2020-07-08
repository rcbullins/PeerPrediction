%Peer prediction inputs script

% % Purpose: Get .mat file containing informationoo needed for peer
%            prediction code

% % Dependencies:StartUp_PeerPrediction
%                PeerPrediction Code
%                PeerPrediction Data
%                Buzcode

% % Input:

% % Output:1) Raster of Graphs Anatomically organized 
%             Raster in Firing Rate order

% % History: 02/04/20 start-RB

% % To do: 
for iSess = 1:length(sessions)

    pathInfo.RecPath     =  [basepath sessions{iSess}];
     
    cd(pathInfo.RecPath);
    
    selecSession    = sessions{iSess};
    datfileName     = [sessions{iSess} '.dat'];
    
    % Define specific session params
    sessionInfo                 = bz_getSessionInfo(cd);
    params.nChans               = sessionInfo.nChannels;
    params.sampFreq             = sessionInfo.rates.wideband;
    params.juxtaspikes_noTM     = 1;
    winRange = 0:5:50;
    
    disp(['Currently evaluating session:' sessions{iSess}])
%% Load in Spikes and lfp
    spikes = bz_GetSpikes;
    %spikes = bz_GetSpikes('sortingMethod','kilosort');
    channel1_lfp = bz_GetLFP(0); % for time of recording reference
    params.length_recording = max(channel1_lfp.timestamps);
%% Anatomical Raster
   % In Channel Order on Probe (given by xml file)
   xlim_lower = 100; % lower limit x (seconds)
   xlim_upper = 101; % upper limit x (seconds)
   PP_rastor_anatom(params, spikes, xlim_lower, xlim_upper)
%% Firing Rate Rastor
   PP_rastor_firingrt(params, spikes, xlim_lower, xlim_upper)
   
%% Convert Timestamps into 1ms bins and assign bins 0/1
   [binned_spikes] = spiketimes2binary(params, spikes.times);

%% ExtraPredictor 1: Peer Prediction + Velocity
   
    cd(pathInfo.RecPath);
   %First get info of analogin
    %read_Intan_RHD2000_file  % click on file
    rhdfilename = ['info.rhd'];
    read_Intan_RHD2000_file_noprompt(rhdfilename)

    analogin_file   = ['analogin.dat'];

    % Get analogin values
    if contains(selecSession,'m1_181220_151127')
        parameters.analoginCh.pulse = 1;
    else
        parameters.analoginCh.pulse = 7;
    end
    parameters.analoginCh.wheel = 2;
    parameters.analoginCh.reward = 1;

    options.downsampleFactor = 30; %sampling to .001ms

    %Get values from _analogin.dat (pulse, water etc.)
    [analogin.pulse, analogin.pos, analogin.reward, analogin.ts] = getAnaloginVals(selecSession,parameters,board_adc_channels,options);
    
    velocities = [diff(analogin.pos)]; % 1 x numBins
    
%% ExtraPreictor 2: Peer Prediction + Position    %Run Velocity Code before this
    position_coords = analogin.pos;
    position_coords(1) = []; %delete first position value
  

%% ExtraPredictor 3:Peer Prediction + Theta Phase %%SAMPLING PROBLEMMMMM******
 cd(pathInfo.RecPath);
  highpass_freq = 6; % passband high 
   lowpass_freq = 10; % passband low 
  theta_channel = 1; %What channel to measure theta on
 % channel_lfp = bz_GetLFP(theta_channel); %use binary 
  channel_lfp = bz_LoadBinary('amplifier.dat','frequency',30000,'nChannels', 64, 'channels',1, 'downsample', 30);
  half_freq = params.sampFreq/2;
  
  % High pass Butter Filter
  [num_coeff_high, denom_coeff_high] = butter(1,highpass_freq/half_freq, 'high');
  % Filter lfp_data through highpass butter filter
  filtered_lfp_high = filtfilt(num_coeff_high, denom_coeff_high, double(channel_lfp));
  % Low Pass Butter Filer
  [num_coeff_low, denom_coeff_low] = butter(1,lowpass_freq/half_freq, 'low');
  %Filter high pass filtered lfp_data through lowpass butter filter
  filtered_theta = filtfilt(num_coeff_low, denom_coeff_low, filtered_lfp_high);
  
  fft_theta = fft(filtered_theta);
  amplitude_theta = abs(fft_theta);
  phase_rad = angle(fft_theta);
  phase_rad(1) = [];
  phase_rad = phase_rad.';
  %convert radians to degrees
  phase_degrees = rad2deg(phase_rad);
 
%% Make mat file of peer prediction inputs
     extraPredictors = [phase_rad;velocities;position_coords];
     cd(pathInfo.RecPath);
    save([selecSession ' peerPrediction_inputs'], 'binned_spikes', 'winRange','extraPredictors','phase_rad','velocities','position_coords','selecSession', '-nocompression', '-v7.3')
 %% 
 disp('Yay science!');

end
