% StartUp_PeerPrediction
%     --Defines paths, parameters, recording sessions,options--
% By:      Reagan  on  02/04/20

% INPUTS: -Drive on Computer with ... -Code and Data
%                                     -Buzcode
%         -Which JuxtaSorter & which ExtraSorter

%% DEFINE DRIVE & SORTERS - EVERYTIME :)
    drive_code = 'C';     
      % Reagan = E 
      % Lianne = E
    drive_buzcode = 'C';
      % Reagan = E
      % Lianne = C
     
%% Add & Define Generic Pathways
    if strcmp(drive_code, 'D')
        addpath(genpath('D:\Data\PeerPrediction_RSC\'))
        addpath(genpath('D:\Code\Github\'))
        basepath = ([drive_code ':\Data\PeerPrediction_RSC\']);
    end
    if strcmp(drive_code, 'E')
         addpath(genpath('C:\Users\lklaver\Documents\GitHub\buzcode'))
        addpath(genpath('E:\Dropbox\Code\PeerPrediction\'))
         addpath(genpath('E:\Code\Github\buzcode-dev\'))
        basepath = ([drive_code ':\Data\PeerPrediction_RSC\']);
    end

%% Recording Sessions Available 
%Sessions availble 

    sessions = {'m115_191203_152410_2', ...
                'm115_191203_081109'}; % Kaiser's data
       
%% Define Generic Params (intrinsic properties of recording)
% index top to bottom anatomically
    params.Probe0idx  = [58 56 54 51 53 47 63 61 48 60 62 52 ...
                         50 55 57 59 32 34 36 38 40 42 44 46 ...
                         49 45 43 41 39 37 35 33 30 28 26 24 ...
                         22 20 18 16 15 19 21 23 25 27 29 31 ...
                         4 6 8 13 11 17 1 3 14 2 0 10 12 9 7 5];
    params.Probeflip    = flip(params.Probe0idx); %0idx 
    params.Probeflip(1) = []; % rm juxta
    params.juxtachan    = 1;
    params.Probe        = params.Probe0idx +1; %1-idx
    params.sampFreq     = 30000;
    params.nChans       = 64;
   
    
%% Define Generic Options   
    %for determining matches
    opts.rangeSpkBin        = .0001; %binsize for extra occurring before or after juxta % maybe a bit wider for james' irc2 output?
    
    %for spectrograms
    opts.timWinWavespec     = 250; %ms
    opts.freqRange          = [1 500];
    opts.numFreqs           = 100;%ops.freqRange(end)-ops.freqRange(1);
    opts.bltimvec           = 10*501-1+250; %5259
        
    %saving and plotting
    opts.SampFreq           = 30000; % should be parameter
    opts.doSave             = 1;
    opts.doPlots            = 0;
    opts.doSaveMat          = 1;
    opts.doSaveSpikes       = 1;
    
    % variable options JuxtaSorter
    opts.intervals          = [0 Inf]; %in sec - change to desired time (find via neuroscope) multiple intervals can be assigned to multiple rows
    opts.downsamplefactor   = 1;
    opts.tempThr            = 15;
    opts.SNRthr             = 7; % figure this one out per cell PARAM SEARCH
    opts.filter             = 'butterworth';
    opts.hpfreq             = 300;
    opts.buttorder          = 1;
    opts.firorder           = 256;
    opts.templateMatch      = 0;
    opts.spikeSamps         = [-40:55];
    
    % loading files
    opts.juxtaFileName = 'juxtaspikes_noTM.mat';
        
    % for CCG    
    opts.ccgBinSize = 0.0015;
    opts.ccgDur     = 0.1;


