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
    [vel_cm_s, time, dt] = getVelocity(analogin, params);
    %gives velocity in 30,000 samples per second smoothed 
    
 %% Find Epochs of running vs not running
    [runEpochs, runIdx] = getRunEpochs(vel_cm_s, dt, time)
    
    total_time_run = 0
    for i = 1:159 %hard coded
        int_row = runEpochs(i,2) - runEpochs(i,1);
        total_time_run = total_time_run + int_row
    end
    total_time_min = total_time_run /60; %28min running out of baseline 40 min
    
    