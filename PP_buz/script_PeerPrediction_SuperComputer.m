%Script to run peer prediction function with
%    input mat file containing: binned_spikes, winRange, extrapredictors,
%           phase_rad, velocities, position_coords, selecSession
%% Adding Paths 
%buzcode & script
    addpath(genpath('/home/reagan4/Code'))
    addpath(genpath('/home/reagan4/Results'))
    %addpath(genpath('/groups/inhibitorange/peerPrediction'))  
          % where data will may need two // in front of groups
%% Load mat file
cd('/home/reagan4/Code') % CHANGE eventually to data folder
load([selecSession ' peerPrediction_inputs.mat')

%% All Together
    cd('/home/reagan4/Code') %necessary?? need to be in path of script?
   [dev devControl] = bz_peerPrediction(binned_spikes(:,1,:), winRange, extraPredictors); %The second trial is just a place holder... only need 1 'trial'   
    cd('/home/reagan4/Results')
    save('ws_full_dev_m115_2', 'params', 'ops', 'dev', 'devControl', '-nocompression', '-v7.3')
   %% Extra Predictors by themselves
   
       cd('/home/reagan4/Code') 
   [dev_phase devControl_phase] = bz_peerPrediction(binned_spikes(:,1,:), winRange, phase_rad);  
       cd('/home/reagan4/Results')   
       save('ws_phase_dev_outputs_m115_2','params', 'ops', 'devControl_phase', 'dev_phase', '-nocompression', '-v7.3')
       
       cd('/home/reagan4/Code') 
   [dev_velocity devControl_velocity] = bz_peerPrediction(binned_spikes(:,1,:), winRange, velocities);  
       cd('/home/reagan4/Results')
       save('ws_velocity_dev_outputs_m115_2', 'params', 'ops', 'devControl_velocity', 'dev_velocity', '-nocompression', '-v7.3')
    
       cd('/home/reagan4/Code') 
   [dev_position devControl_position] = bz_peerPrediction(binned_spikes(:,1,:), winRange, position_coords);  
       cd('/home/reagan4/Results')   
       save('ws_position_dev_outputs_m115_2', 'params', 'ops', 'devControl_position', 'dev_position', '-nocompression', '-v7.3')
 