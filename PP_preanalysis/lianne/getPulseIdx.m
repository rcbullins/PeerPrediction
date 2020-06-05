function [pulseIdx, noPulseIdx, pulseEpochs] = getPulseIdx(pulse, analogin, tr_ep)

% now depends on trials, make it more robust for all pulses

ts = analogin.ts;

%% Finding depolarizations over time

diffPulse = diff(pulse);
posPulseIdx = diffPulse > 0.2; % now hardcoded
negPulseIdx = diffPulse < -.2; % now hardcoded

selPosPulseIdx  = find(posPulseIdx~=0);
selNegPulseIdx  = find(negPulseIdx~=0);
pulseEpochs     = [ts(selPosPulseIdx+1)' ts(selNegPulseIdx+1)']; % +1 , because diff % pulse Time

% tr_ep; % contains the start and stop of trial
% trlIdx of depolarization falls within trial y/n


countDep = 0;
countNotDep = 0;

for iPulse = 1:length(pulseEpochs)
    countDep = countDep+1;
    for iTr = 1:length(tr_ep)
        if pulseEpochs(iPulse,1) > tr_ep(iTr,1) && pulseEpochs(iPulse,2) < tr_ep(iTr,2)
            pulseIdx(countDep) = iTr;
            break
        end
    end
end

trialVec = 1:length(tr_ep);
noPulseIdx = trialVec(~ismember(trialVec,pulseIdx));

end
