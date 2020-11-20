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
    data_path = [basepath 'PP_RSC_Data\' session_name];
%Add Paths
    addpath(genpath(data_path));
    addpath(genpath([basepath 'Code\']));
    addpath(genpath([basepath 'Sam_Code\']));
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
    min_thresh = .75
    [runEpochs, runIdx] = getRunEpochs(vel_cm_s, dt, time, min_thresh)
    
     %find how long each epoch last
     length_run = zeros(length(runEpochs),1);
     for i = 1:length(runEpochs)
         length_run(i,1) = runEpochs(i,2)-runEpochs(i,1);
     end
     % find epochs greater than 5 seconds
     time_thr = 5;
     long_run_epochs_idx = find(length_run(:,1) >= time_thr)
     runEpochs_long = runEpochs(long_run_epochs_idx,:);
     % sanity check - correct epochs
     figure
     plot(time(2:end),vel_cm_s);
     xlim([0 4500]);
     title('vel_cm_s');
     hold on
    %%
    total_time_run = 0
    for i = 1:159 %hard coded
        int_row = runEpochs(i,2) - runEpochs(i,1);
        total_time_run = total_time_run + int_row
    end
    total_time_min = total_time_run /60; %28min running out of baseline 40 min
    
    length_run = zeros(1120,1);
    for i = 1:1120
         length_run(i,1) = runEpochs(i,2)-runEpochs(i,1);
    end
    
 %% Find Epochs of not running

 noRun = zeros(455,2)
 for i = 1:454
     noRun(i,1) = runEpochs(i,2)+.0001
     noRun(i,2) = runEpochs(i+1,1)-.0001
     
 end
 
 NR_epochLength = zeros(455,1)
 
 for i = 1:159
     NR_epochLength(i,1) = noRun(i,2)-noRun(i,1)
     
 end
 %% Firing Rate
 
 new_spikes = {};
        for i = 1:19
            new_spikes{i} = spikes.times{supergoodUnitQualIdx(i)}
        end
        
        firing_rate = zeros(1,19)
      for i = 1
        length_time = new_spikes{i}(length(new_spikes{i})) - new_spikes{i}(1);
        num_spikes = length(new_spikes{i});
        firing_rate(1,i) = num_spikes/length_time;
      end
 %% Get power spectrum
 channel = 57;
 %epoch = [1547.5 2143.4];
 epoch = [1208 1309] %136
 cd(data_path)
 getPowerSpectrum(channel, epoch, session_name)
 
 