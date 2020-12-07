%% PP_assembly_specifications

%Add Basepath for all code and data used
    basepath = ('C:\Users\rcbul\Documents\English Lab\');
%Define Recording Session Name
    %session_name = 'm115_191203_152410_n';
    session_name = 'u19_200313_155505'; %NUM 1 -running
    %session_name = 'u19_200310_135409'; %NUM 2
    %session_name = 'u19_200313_120452'; %NUM 3
    %session_name = 'u21_200305_153604'; %NUM 4 -runninng
    %session_name = 'u21_200309_142534'; %NUM 5
    %session_name = 'u26_200306_172032'; %NUM 6
%Deine DataPath that contains list of session names;
    data_path = [basepath 'PP_Data\' session_name];
%Add Paths
    addpath(genpath(data_path));
    addpath(genpath([basepath 'Code\']));
    addpath(genpath([basepath 'buzcode-dev\']));
%% Set Graph Defaults Now
    SetGraphDefaults;
%% Params
    params.radiusDisk   = 26;
    params.circDisk     = 2*pi*params.radiusDisk;
%%  Get Velocity
    cd(data_path)
    load([session_name '_analogin.mat']);
    [vel_cm_s, time, dt] = getVelocity(analogin, params, session_name);
    %gives velocity in 30,000 samples per second smoothed 
    
 %% Find Epochs of running 
 %make .75
    min_thresh = 5
    [runEpochs, runIdx] = getRunEpochs(vel_cm_s, dt, time, min_thresh)
    
     %find how long each epoch last
     length_run = zeros(length(runEpochs),1);
     for i = 1:length(runEpochs)
         length_run(i,1) = runEpochs(i,2)-runEpochs(i,1);
     end
     % find epochs greater than 5 seconds
     time_thr = 5;
     runIdx_long = find(length_run(:,1) >= time_thr)
     runEpochs_long = runEpochs(runIdx_long,:);
     % sanity check - correct epochs
     figure
     plot(time(2:end),vel_cm_s);
     xlim([0 4500]);
     title('vel_cm_s');
     hold on
 %% Find Epochs of not running

 noRunEpochs = zeros(length(runEpochs),2)
 for i = 1:length(runEpochs)-1
     noRunEpochs(i,1) = runEpochs(i,2)+.0001
     noRunEpochs(i,2) = runEpochs(i+1,1)-.0001
 end
 
 length_noRun = zeros(length(noRunEpochs),1)
 
 for i = 1:length(noRunEpochs)-1
     length_noRun(i,1) = noRunEpochs(i,2)-noRunEpochs(i,1)
 end
     time_thr = 5;
     noRunIdx_long = find(length_noRun(:,1) >= time_thr)
     noRunEpochs_long = noRunEpochs(noRunIdx_long,:);
 
 %% Get Power Spectrum
 rippleChan = 57;
 epoch = runEpochs_long(2,:) %136
 cd(data_path)
 getPowerSpectrum(rippleChan,epoch)

 %% Alt Get wavespec (have mat file of data structure)
%  runWaveSpec(opts.rippleChan, runEpochs_long, runIdx_long)

   % INPUT: Define which epoch you want to graph
    epoch_idx = 4; %any epoch in runEpochs_long
 
    %Define the start and stop of the epoch you want to graph
        epoch_lim = runEpochs_long(epoch_idx,:); %seconds    
    %Plot Velocity 
        figure
        subplot(2,1,1)
        t = time;
        t_ds = downsample(t,10);
        vel_ds = downsample(vel_cm_s, 10);
        plot(t_ds, vel_ds)
        xlim([epoch_lim(1)-10 epoch_lim(2)+10])
        xlabel('Time (s)')
        ylabel('Velocity')
        hold on
        xline(epoch_lim(1))
        xline(epoch_lim(2))
    
    % Wavespec plot
        cd(data_path)
        load('wavespecall.mat'); %sampled at 1250 per second
        subplot(2,1,2)
        normdata = abs(ws_temp.data);
        %normdata_ds = downsample(normdata, 125); %10 samples per second
        imagesc(normdata');
        set(gca,'YDir','normal');
        colormap(jet)
        xlabel('Time(s)');
        ylabel('Frequency(Hz)');
        hold on
        xline(epoch_lim(1)*1250)
        xline(epoch_lim(2)*1250)
        base_time = 10*1250;
        xlim([epoch_lim(1)*1250-base_time epoch_lim(2)*1250+base_time])
 %% Wavespec average plot for all run epochs
   %how many seconds to plot before and after run
   length_t = 7; %seconds
   start_time = runEpochs_long(:,1) - length_t;
   end_time = runEpochs_long(:,1) + length_t;
   start_sample = start_time*1250;
   end_sample = end_time*1250;
   
   cd(data_path)
   load('wavespecall.mat'); %sampled at 1250 per second  
   normdata = abs(ws_temp.data);
   ws_data = {};
   ws_avg = {};
   for ii = 1:length(runEpochs_long)
       ws_data{ii} = normdata(start_sample(ii):end_sample(ii)-1,:)
   end
   %average by element over array for each frequency and time point
   ws_cat = cat(3,ws_data{:});
   ws_avg = mean(ws_cat,3);
   
   imagesc(ws_avg');
   set(gca,'YDir','normal');
   colormap(jet)
   xlabel('Time(s)');
   ylabel('Frequency(Hz)');
   xticks(1:8749:17500)
   xticklabels({'-7','0','7'})
  