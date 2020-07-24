%Script for Analyzing PP Result Data

%Inputs: Peer Prediction function inputs and outputs

%        See first section of code to define...
%           
%        Basepath: Containing all Code and Data Necessary
%        Session Name: Recording Session 
%        DataPath: Path to defined Session Name
%        Result Data Path: Path to outputs from bz_peerPrediction (dev and
%        devControl)
%        Data Format: bz_peerPrediction using link function: Normal or
%        Poisson

%Dependencies: Peer Prediction Function Input & Output Data
%              Raw Spike Struct Data

%Outputs: Assembly Strength
%         Deviance Graphs per cell pair
%         Histogram of Optimal PP Window
%         ScatterPlot of Error vs Optimal Time Window

% Created: 3/23/20 by Reagan Bullins
% Updated: 6/01/20 by Reagan Bullins  
% Updated: 6/29/20 by Reagan Bullins

% NOTE: Change paths to adapt to your computer
%       Go to 'EXTRA Specifications' Section to specify ... 
%           subset of data to analyze (or all)
%           pairs you want to graph (or all- but that's a lot)
%           winRange

%% Temporary for Test running code
winRange = [1 2 4 8 16 32 64 128 256 512 1024]; %add 0
%winRange = [1:150];
tic
[dev devControl] = bz_peerPredictionRB(binned_spikes(:,1,:),winRange,[],pairsToRun(300,:));
disp('Done')
cd('C:\Users\rcbul\Documents\English Lab\PP_RSC_Data\Testing\WholeRec');
save('pair300_bin1_shuff3.mat','dev', 'devControl')
toc
%% Define Paths: 
% FORMAT: Data should be in folders by session name.
%         Each Session Folder will have...
%            PeerPrediction results folder (dev and devControl)
%            pairsToRun
%            spikes.cellinfo.mat
%            binned_spikes

%Add Basepath for all code and data used
    basepath = ('C:\Users\rcbul\Documents\English Lab\');
    addpath(genpath(basepath));
%Define Recording Session Name
    %session_name = 'm115_191203_152410_n';
    session_name = 'u19_200313_155505';
%Deine DataPath that contains list of session names;
    data_path = [basepath 'PP_RSC_Data\' session_name];
%Define ResultPath that contains results from peer prediction function
    %data_set = 'PeerPrediction_binLog';
    %data_set = 'PeerPrediction_bin2';
    data_set = 'Buzcode'
    result_data_path = [data_path '\' data_set '\'];
%Result data format:NORMAL OR POISSON
    %data_format = 'pp_batch' ;
     data_format = 'pp_poisson'; %FOR NOW ALL ARE Poissons
    
%% EXTRA Specifications
  % Analyze all pairs or a subset?
        choice_analysisPairs = "All";
         %choice_analysisPairs = [1:25];
        %choice_analysisPairs = "WeakPairs";
        %choice_analysisPairs = "StrongPairs";
  % Graping (Deviance, Cross Corr, Rastor): Which pairs?
        choice_graphPairs = [1];
  % What winRange was used?
        %winRange = (0:3:150);
        %winRange = (0:150);
        winRange = [1 2:2:30 32:5:64 128 256]; %used for pp
        %winRange = [1:2:64]; %bin2
         %winRange = [1 2:2:10 12:1:64 66:4:82]; %mixLog
        %winRange = [1:1:70]; %1to70
  
%% Graphing Defaults
SetGraphDefaults;

%% Load and Concatenate Data

[dev, devControl] = PP_Load_and_Concatenate(result_data_path);

%% If Win Range is NOT continuous 
% if the third index in winRange is not the second index in winRange + 1,
% then the winRange must not be continuous. Therefore take out all rows
% with all zeros
if winRange(3) ~= winRange(2)+1
    %take away all zero rows in dev and devControl if window isn't
    %sequential
     dev = dev(any(dev,2),:);
     winRangeT = winRange +1;
     devControl = devControl(winRangeT,:,:);
end

%% Assembly StrengthOnly 

[dev_smoothed, ratio_strength, dev_min, weak_pairs, strong_pairs] = PP_AssemblyStrength (dev, devControl, data_format, winRange);

%% Deviance Analysis 

[pairs_for_analysis, min_win_total, min_win_pairs] = PP_DevianceAnalysis (dev_smoothed ,dev_min, choice_analysisPairs, weak_pairs, strong_pairs, winRange);

%% Deviance Graph 

PP_DevianceGraphs(dev, dev_smoothed, devControl, dev_min, ratio_strength, choice_graphPairs, winRange);

%% Histogram of Optimal Time Windows : ALL
bin_win_count = 1; %binned within 3 ms time windows -- CAN adjust
bin_win_max = winRange(length(winRange));
bin_win_max = bin_win_max + bin_win_count;

[optimal_window] = PP_TimeWindow_Histogram (bin_win_count, bin_win_max, min_win_total); %optimal window is median

%% Histogram of Optimal Time Windows: Analysis Specifications
% bin_win_count = 1; %binned within 3 ms time windows -- CAN adjust
% bin_win_max = winRange(length(winRange)); % CAN adjust 
% bin_win_max = bin_win_max + bin_win_count;

[optimal_window] = PP_TimeWindow_Histogram (winRange, min_win_pairs'); %optimal window is median

%% Histogram of Assembly Strength

bin_interval = 2; 
PP_AssemblyStrength_Hist (ratio_strength, bin_interval);

%% Scatterplot of Optimal Time Windows: All

y_plot_limit = 20;
PP_OptimalWindow_ScatterPlot(min_win_total,dev_smoothed,devControl, y_plot_limit, winRange);

%% Scatterplot of Optimal Time Windows: Analysis Specifications

y_plot_limit = 20;
PP_OptimalWindow_ScatterPlot(min_win_pairs,dev_smoothed,devControl,y_plot_limit, winRange);

%% Cross Correlograms of Pairs :Graph specifications
  % Load in data 
     cd(data_path)
     load('pairsToRun.mat')
     load('binned_spikes.mat')
     load([session_name '.spikes.cellinfo.mat']);
  % Get length of recording
    [~, ~, bin_win_max] = size(binned_spikes);
    clear binned_spikes
  % Define bin size
    bin_size = .01; %default is .01
    PP_crosscorr(bin_size, spikes, choice_graphPairs, pairsToRun, bin_win_max)

%% Rastor for Actual and Predictor: Graph specifications

  % Load in Data
    cd(data_path)
    load('pairsToRun.mat')
    load('binned_spikes.mat')
    load([session_name '.spikes.cellinfo.mat']);
  % Define second you want to graph and pair
    choice_sec = 105; %which second you want to plot
    graph_pair = 1; %idk why only one graph at a time works.
    win_num = [24 52 128 256]; %the ms windows we will look at FOR SPEED

[smoothedTrains] = PP_Raster_SmoothedTrains(graph_pair, choice_sec, binned_spikes, pairsToRun, spikes, win_num)

%% Consistent Spiking Visualize 
% in sam's code folder
cd(data_path)

[spikeRates] = PP_dotplot_spikerate(binned_spikes);

%% Place Field 


%% Code for seeing cell specific involvement per bin size
list_BIG_pairs = zeros(length(list_BIG),2);
for i = 1:length(list_BIG)
    list_BIG_pairs(i,1) = pairsToRun(list_BIG(i),1);
    list_BIG_pairs(i,2) = pairsToRun(list_BIG(i),2);
end

list_big_pairs = zeros(length(list_big),2);
for i = 1:length(list_big)
   list_big_pairs(i,1) = pairsToRun(list_big(i),1);
   list_big_pairs(i,2) = pairsToRun(list_big(i),2);
end

hist_BIG = zeros(580,1)
hist_BIG(1:290,1) = list_BIG_pairs(:,1);
hist_BIG(291:580,1) = list_BIG_pairs(:,2);

hist_big = zeros(150,1)
hist_big(1:75,1) = list_big_pairs(:,1);
hist_big(76:150,1) = list_big_pairs(:,2);

histogram(hist_big)

