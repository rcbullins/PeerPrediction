
%% normalize the waveform between 1 and -1
for iUnit = 1:length(spikes.UID)
    if isempty(spikes.rawWaveform{iUnit})
        continue
    else
        
        nSamp = size(spikes.rawWaveform{iUnit},1);
        spkTrace = spikes.rawWaveform{iUnit};
        [minWave, minIndx] = min(spkTrace);
        
        [maxWave, maxIndx] = max(spkTrace(minIndx:end));
        maxIndx = maxIndx + minIndx;
        
        normWave = 2*((spkTrace - minWave)./ (maxWave-minWave))-1;
        
        
        %% shift all peaks to the same sample
        
%         peakSamp = 61; %  61th sampls is the peak value
%         shiftIndx = minIndx-peakSamp;
%         minShift = min(shiftIndx);
%         
%         %% vv Hier ben je mee bezig, Lianne
%         padNormWave = [nan(size(normWave,1),abs(minShift)) normWave nan(size(normWave,1),abs(maxShift))];
%         
%         nUnits = size(padNormWave,1);
%         shiftWave = nan(nUnits,nSamp);
%         for iUnit = 1:nUnits
%           shiftWave(iUnit,:) = padNormWave(iUnit,(maxIndx(iUnit)-(peakSamp-1))+abs(minShift): ...
%             (maxIndx(iUnit)+24)+abs(minShift));
%         end
%         
%         figure; plot(shiftWave')
        %% ^^ Dit werkt nog niet optimaal
        
        shiftWave(iUnit,:) = normWave;

    end

%% Peak to through dist (waveform dur in Barrel paper). Seems like Martin interpolated the waveforms
TroughtoPeakDur(iUnit) = (maxIndx - minIndx)*1/nSamp; %30 samples
end
figure,plot(shiftWave')
 
figure; hist(abs(TroughtoPeakDur),[0:1/30:2.5]) % 30 samples per wave

figure,plot(shiftWave','b')
hold on
plot(shiftWave(TroughtoPeakDur<0.2,:)','k') % 0.1 for RSC

INTIndx= TroughtoPeakDur<0.1


%% Repolarization at .45ms (2nd measure from the barrel paper)
% samp 8 = 0time
% .45ms is 14.4 samples.. which is a little strange, so I guess take sample 8+14 = 22
% .55ms is 17.6 samples 8+18 = 26

% Lianne: for my data
% 121/30000 = 4ms window in total

repolarization = shiftWave(:,20); % what's this?
figure; hist(repolarization,50)

ind1 = repolarization <= -.85;

%% 
figure; hold on;
plot(shiftWave(ind1,:)','k')


% % %% FROM HERE STILL TO WRITE
% % %% Some ugly units exclusion
% % 
% % figure; plot(repolarization, TroughtoPeakDur,'o')
% % figure; plot(log10(freq.rec), TroughtoPeakDur,'o')
% % 
% % % %% Find some of the ugly shapes
% % % 
% % % bInd1 = (shiftWave(:,3) < -0.9);
% % % bInd2 = (shiftWave(:,14) > 0.6);
% % % bInd3 = (shiftWave(:,19) > -0.65 & shiftWave(:,19) < -0.4);
% % % bInd4 = (shiftWave(:,15) > -0.37 & shiftWave(:,15) < -0.15);
% bInd5 = (shiftWave(:,26) < -0.98);
% 
% badInd = bInd1 | bInd2 | bInd3 | bInd4 | bInd5;
% 
% figure; hold on
% plot(shiftWave(~badInd,:)','k')
% plot(shiftWave(badInd,:)','r')
% plot(shiftWave(70,:)','g')
% 
% find(badInd)


%% Some firing rate inclusion