%% Define Varibables
    VT_PID = ''; %virginia tech PID
    session_name = ''; %session name you want to run

%% Add Globus Paths 
    addpath(genpath(['/home/' VT_PID '/Code']))
    addpath(genpath(['/home/' VT_PID '/Results']))
    addpath(genpath(['/home/' VT_PID '/Data']))

%% Load input mat file
    cd(['/home/' VT_PID '/Data'])
    load([session_name '_globus.mat']) %need spikes struct, bin_list, session_name, results folder,epochRun
%% Run Script
cd(['/home/' VT_PID '/Results/' session_name '/' results_folder]);
%run sixth set of 6
tic
for ibin = 28:33  %length(bin_list)
    [log_likelihood,weights] = CrossValidationAssemblyPrediction(spikes, 'dt', bin_list(ibin), 'epoch', epochRun);
    save([num2str(bin_list(ibin)) '_bn_asmb.mat'], 'log_likelihood', 'weights', '-v7.3');
    clear log_likelihood weights
end
toc
