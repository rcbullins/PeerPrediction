addpath(genpath('E:\Dropbox\Code\english_lab\'))
addpath(genpath('E:\Dropbox\Code\buzcode\'))
addpath('E:\Dropbox\Code\intan')
addpath('E:\Dropbox\Code\buzsupport')
%%
basename        = 'm1_181220_151127';  % m121_191126_124321
basepath        = 'E:\Data\interneuron_stim\';
% 
% basename ='m121_191126_124321'
% basepath = 'E:\Data\odor_pilot\m121\'

% basename = 'm122_191202_144907'
% basepath = 'E:\Data\odor_pilot\m122\'

analogin_path   = fullfile([basepath filesep basename]);
spike_path      = fullfile([basepath filesep basename filesep 'ks2ec']);

cd(spike_path)

%Load the spikes
spikes = bz_LoadPhy; %only works when clusters are labeled 'good'

cd(analogin_path)

%% Split Data into Trials
% First get info of analogin
% read_Intan_RHD2000_file  % click on file
rhdfilename = [basename '_info.rhd'];
read_Intan_RHD2000_file_noprompt(rhdfilename)

analogin_file   = [basename, '_analogin.dat'];

% Get analogin values
if contains(basename,'m1_181220_151127')
    parameters.analoginCh.pulse = 1; 
else
    parameters.analoginCh.pulse = 7;
end
parameters.analoginCh.wheel = 2;
parameters.analoginCh.reward = 1;

options.downsampleFactor = 300; %(samplingRate/100Hz to get 100Hz)

%Get values from _analogin.dat (pulse, water etc.)
[analogin.pulse, analogin.pos, analogin.reward, analogin.ts] = getAnaloginVals(basename,parameters,board_adc_channels,options);

% Separate the wheel in trials
[len_ep, ts_ep, vel_ep, tr_ep, len_ep_fast, ts_ep_fast, vel_ep_fast] = getWheelTrials(analogin);


% % % pulseEpochs = [ts(selPosPulseIdx+1)' ts(selNegPulseIdx+1)']
%%
for iUnit = 1:length(spikes.UID)
    
   [status,interval] = InIntervals(spikes.times{iUnit},tr_ep);
    
    
    % spk_ep %30000 Hz
    % find out what position each spike is
    % by binning positions
    
    trialTimes  = ts_ep; % 100 Hz
    trialPos    = len_ep; % 100Hz
    
    
    
    for iTr = 1:length(len_ep)
        
%         spk_ep{iUnit}.trial{iTr}  = spikes.times{iUnit}(interval==iTr);
        
        lengthWheel = 1.2; % in velocity measure %double check this ltaer % 1.2 for m1
        posBinsCount = histoc(len_ep{iTr,1},0:0.05:lengthWheel); % 6 position bins to test 
        
        ts_posbin = mat2cell(trialTimes{iTr},posBinsCount);% all timestamps per trial for 6 positions.

         for iPosBin = 1:length(posBinsCount) % hardcoded now
            if ~isempty(ts_posbin{iPosBin})
                posbinTs(iPosBin,:) = [ts_posbin{iPosBin}(1) ts_posbin{iPosBin}(end)]; % in time
                % seconds in bin per trial
                secTrPosBin = [posbinTs(iPosBin,2) - posbinTs(iPosBin,1)]; % in sec
            else
                posbinTs(iPosBin,:) = [0.0000001 0.0000002];
                secTrPosBin = [0.0000000001];
            end
        end
        
        
        %how many spikes in these intervals
        [status, interval]= InIntervals(spikes.times{iUnit},posbinTs);
        
        % now find out how many spikes/s are within these posbinS. (normalized FR 0 to 1) By
        % calculating length posbinTs (end-1) and num timestamps.
        numSpkBin      = histoc(interval(interval>0),1:size(posbinTs,1)); % how many timestamps or position values should there be per bin
        
        
        spkPerSec{iUnit}.trial{iTr} = numSpkBin/secTrPosBin;
    end
    % time stamps per bin
    
    unitRM_Matrix = cell2mat(spkPerSec{iUnit}.trial)';
    rowmin = min(unitRM_Matrix,[],2);
    rowmax = max(unitRM_Matrix, [],2);
    resc_unitRM{iUnit} = rescale(unitRM_Matrix,'InputMin',rowmin, 'InputMax', rowmax); % normalized FRs
    
    % figure, imagesc(resc_unitRM{iUnit});
    % clear unitRM_Matrix  rowmin rowmax numSpkBin posbinTs
end

% this does not work now, obtained pulseEpochs are not correct

[depIdx, nonDepIdx, pulseEpochs] = getPulseIdx(analogin.pulse, analogin, tr_ep);%parameters.analoginCh


%% Plot ratemaps per unit (each row is 1 trial) Depolarizations vs non-depolarization trials

for iUnit=1:length(spikes.UID)
    figure
    subplot(2,1,1)
    imagesc(resc_unitRM{iUnit}(depIdx,:))
    box off
    set(gca,'TickDir','out')
    xlabel('position (bins)')
    ylabel('trials')
    c1 = colorbar;
    c1.Label.String = 'normalized FR';
    title('Normalized FR for trials with depolarization')
    
    subplot(2,1,2)
    imagesc(resc_unitRM{iUnit}(nonDepIdx,:))
    box off
    set(gca,'TickDir','out')
    xlabel('position (bins)')
    ylabel('trials')
    c2 = colorbar;
    c2.Label.String = 'normalized FR';
    title('Normalized FR for trials without depolarization')
end


%% Where is the depolarization? (Location)
% find corresponding location
% pulseEpochs(:,1) %start times of pulse
% depIdx % matching indices of those pulse times
% len_ep % position per trial
% ts_ep % corresponding timestamp
% for iTr=1:length(tr_ep),a(iTr)=sum(ismember(ts_ep{iTr}, pulseEpochs(:,1))), end

%% Rasters per unit Trials with Depolarizations, Without Depolarization
% deze doet het niet bepaald nu. 

for iUnit = 1:length(spikes.UID)
    % Where are the spike times per trial?
    [status, interval]= InIntervals(spikes.times{iUnit},tr_ep);
    numSpkTrial = histoc(interval(interval>0),1:size(tr_ep,1));
    spike_ep{iUnit}.trial  = mat2cell(spikes.times{iUnit}(status),numSpkTrial); % spiketimes per trial
end
%
%plot rasters
for iUnit = 1:length(spikes.UID)
    figure
    plotSpkOffset = 0;
    for iTr = 1:length(tr_ep)
        %     for iSpk = 1:length(spike_ep{iUnit}.trial{iTr})
        plotSpikeTs = spike_ep{iUnit}.trial{iTr} - min(spike_ep{iUnit}.trial{iTr}); % to align all trials subtract the first ts
        if isempty(plotSpikeTs)
            break
        else
            line([plotSpikeTs plotSpikeTs],[plotSpkOffset plotSpkOffset+1],'Color','red'),
            plotSpkOffset = plotSpkOffset-1;
            hold on
        end
    end
end
%% These are also rasters and these do work 
for iUnit = 1:length(spikes.UID)
    figure
    plotSpkOffset = 0;
    countDepRas =0;
    subplot(2,1,1)
    for iTr = depIdx
        countDepRas = countDepRas + 1;

        plotSpikeTs = spike_ep{iUnit}.trial{iTr} - min(spike_ep{iUnit}.trial{iTr}); % to align all trials subtract the first ts
        if isempty(plotSpikeTs)
            break
        else
            line([plotSpikeTs plotSpikeTs],[plotSpkOffset plotSpkOffset+1],'Color','black'),
            line([pulseEpochs(countDepRas,1)-min(spike_ep{iUnit}.trial{iTr}) pulseEpochs(countDepRas)-min(spike_ep{iUnit}.trial{iTr})],[plotSpkOffset plotSpkOffset+1],'Color','red')
            plotSpkOffset = plotSpkOffset-1;
            hold on
        end
        xlabel('time')
        ylabel('trials')
        title('rasterplots per trial, red = depol time')
        box off
        set(gca,'TickDir','out')
    end
    
    subplot(2,1,2)
    for iTr = nonDepIdx
           plotSpkOffset = 0;

        plotSpikeTs = spike_ep{iUnit}.trial{iTr} - min(spike_ep{iUnit}.trial{iTr}); % to align all trials subtract the first ts
        if isempty(plotSpikeTs)
            break
        else
            line([plotSpikeTs plotSpikeTs],[plotSpkOffset plotSpkOffset+1],'Color','black'),
            plotSpkOffset = plotSpkOffset-1;
            hold on
        end
        xlabel('time')
        ylabel('trials')
        title('rasterplots per trial, red = depol time')
        box off
        set(gca,'TickDir','out')
    end
end

% depIdx
% nonDepIdx
%% Calulate Alignment to Depolarization

timeBefore = 1;
timeAfter = 1;

trlCenteredDepStart = pulseEpochs(:,1)-timeBefore;
trlCenteredDepStop = pulseEpochs(:,1)+timeAfter;

trlCenteredDep = [trlCenteredDepStart trlCenteredDepStop];


for iUnit  = 1:length(spikes.UID)
    [status, interval]= InIntervals(spikes.times{iUnit},trlCenteredDep);
    numSpkCoD = histoc(interval(interval>0),1:size(trlCenteredDep,1));
    spike_ep_CoD{iUnit}.depol  = mat2cell(spikes.times{iUnit}(status),numSpkCoD);
end

% Rasters aligned to depolarization

for iUnit = 1:length(spikes.UID) 
    figure
    plotSpkOffset = 0;
    countDepRas =0;
    
    depTrTs = tr_ep(depIdx,:);
    subtrTrStart = depTrTs(:,1);
    
    %pulseEpochs(:,1) % pulse start Times
    for iPulse = 1:length(depIdx)
        plotSpikeforHisto{iPulse} = spike_ep_CoD{iUnit}.depol{iPulse} - pulseEpochs(iPulse,1);
        
        countDepRas = countDepRas + 1;
        %         plotSpikeTs = spike_ep_CoD{iUnit}.depol{iPulse} -subtrTrStart(iPulse); % min(spike_ep_CoD{iUnit}.depol{iPulse}) to align all trials subtract the first ts
        plotSpikeTs = spike_ep_CoD{iUnit}.depol{iPulse} - pulseEpochs(iPulse,1); %
        
        
        if isempty(plotSpikeTs)
            continue
        else
            subplot(2,1,1)
            hist(cell2mat(plotSpikeforHisto'),50)
            title(['Unit' num2str(iUnit)])
            xlabel('time(s)')
            ylabel('count')
            ymax = get(gca,'YLim');
            line([pulseEpochs(countDepRas,1)-pulseEpochs(iPulse,1) pulseEpochs(countDepRas,1)-pulseEpochs(iPulse,1)],[0 ymax(2)],'Color','red')
            line([pulseEpochs(countDepRas,2)-pulseEpochs(iPulse,1) pulseEpochs(countDepRas,2)-pulseEpochs(iPulse,1)],[0 ymax(2)],'Color','blue')
            box off
            set(gca,'TickDir','out')
            
            
            subplot(2,1,2)
            line([plotSpikeTs plotSpikeTs],[plotSpkOffset plotSpkOffset+1],'Color','black'),
            line([pulseEpochs(countDepRas,1)-pulseEpochs(iPulse,1) pulseEpochs(countDepRas,1)-pulseEpochs(iPulse,1)],[plotSpkOffset plotSpkOffset+1],'Color','red')
            line([pulseEpochs(countDepRas,2)-pulseEpochs(iPulse,1) pulseEpochs(countDepRas,2)-pulseEpochs(iPulse,1)],[plotSpkOffset plotSpkOffset+1],'Color','blue')
            
            plotSpkOffset = plotSpkOffset-1;
            hold on
            xlabel('time (s)')
            ylabel('trials')
            box off
            set(gca,'TickDir','out')
        end 
    end
end

% double check of alignment is now most elegant way of doing this,
% realign trials to start each trial at 0 perhaps needs to happen earlier
%%
for iUnit = 1:length(spikes.UID)
    figure,
    for iPulse = 1:length(depIdx)
        plotSpikeforHisto{iPulse} = spike_ep_CoD{iUnit}.depol{iPulse} - pulseEpochs(iPulse,1);
    end
    
    hist(cell2mat(plotSpikeforHisto'),50)
end


%%
countDep = 0;
countNonDep = 0;

for iUnit = 1:length(spikes.UID)
    countDep = countDep+1;
    allSpkTrls = spike_ep{iUnit}.trial;
    for iPulse = 1:length(depIdx)
        spiketimes_depol{iUnit}.trial{countDep} = allSpkTrls{iPulse};
    end
    
    for iNoPulse=1:length(nonDepIdx)
        countNonDep = countNonDep+1;
        spiketimes_nopulse{iUnit}.trial{countNonDep} = allSpkTrls{iNoPulse};
    end
    
    spiketimes_depol_allTr{iUnit} = cell2mat(spiketimes_depol{iUnit}.trial');
    spiketimes_nopulse_allTr{iUnit} = cell2mat(spiketimes_nopulse{iUnit}.trial');
end


%%
%CCG
params.sampFreq=30000;

plotops.ccgBinSize  = 0.001;
plotops.ccgDur      = 0.2;
%histogram
plotops.histBinSize = 10;

[ccg_nopulse,t_nopulse] = CCG(spiketimes_nopulse_allTr,[],'Fs',params.sampFreq, 'binSize',plotops.ccgBinSize,'duration', plotops.ccgDur,'norm','counts');
[ccg_depol,t_depol] = CCG(spiketimes_depol_allTr,[],'Fs',params.sampFreq, 'binSize',plotops.ccgBinSize,'duration', plotops.ccgDur,'norm','counts');

%%
figure
plotCount = 1;
for idx_hMFR = 1:10%size(ccg_nopulse,2)
    for iPair = 1:10%size(ccg_nopulse,2)
        subplot(size(ccg_nopulse,2),size(ccg_nopulse,2),plotCount)%size(ccg_nopulse,2),
        plotCount = plotCount+1;
        if idx_hMFR == iPair
            bar(t_nopulse,ccg_nopulse(:,idx_hMFR,iPair),'k')
        else
            bar(t_nopulse,ccg_nopulse(:,idx_hMFR,iPair))
        end
    end
end

%
figure
plotCount = 1;
for idx_hMFR = 1:10%size(ccg_depol,2)
    for iPair = 1:10%size(ccg_depol,2)
        subplot(size(ccg_nopulse,2),size(ccg_nopulse,2),plotCount);%size(ccg_depol,2),size(ccg_depol,2)
        plotCount = plotCount+1
        if idx_hMFR == iPair
            bar(t_depol,ccg_depol(:,idx_hMFR,iPair),'k')
        else
            bar(t_depol,ccg_depol(:,idx_hMFR,iPair))
        end
    end
end
