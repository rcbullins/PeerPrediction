function [len_ep, ts_ep, vel_ep, tr_ep, len_ep_fast, ts_ep_fast, vel_ep_fast] =getWheelTrials(analogin)

% Lianne 20191126, adapted from all the permutations of getDiskTrials
% Reagan 20200218 adapt to fix start and stop intervals to be same length
pos = analogin.pos;
ts = analogin.ts;

%% Finding position over time
% This section is trying with findpeaks

k = normpdf([1:40],20,5);
thres = .5e-3;

% Smooth position (pos) for detection of slope, etc
pos = nanconvn(pos,k);

% figure,
% plot(ts,pos)
% 
% figure,
% plot(ts,pos)

[~, trP_idx] = findpeaks(pos, 'MinPeakProminence', 1, 'MinPeakHeight', 1);
[~, trN_idx] = findpeaks(-pos, 'MinPeakProminence', .5, 'MinPeakHeight', -1);

% figure,
% findpeaks(pos, 'MinPeakProminence', 1, 'MinPeakHeight', 1);
% figure,
% findpeaks(-pos, 'MinPeakProminence', .5, 'MinPeakHeight', -1);

    %If the Start Time (P) vector is longer than the Stop (N), the
    %recording ended with an extra start. So remove last index.
        if length(trP_idx) > length(trN_idx)
            trP_idx(end) = [];
        end
    %If the Stop Time (N) vecotr is longer than the Srart(P), the recording
    %began with an extra stop. So remove the first index.
        if length(trP_idx) < length(trN_idx)
            trN_idx(1) = [];
        end
        
% get trials timestamps Positive and Negative (= begin and end of a lap)
% pad_trN_idx = [trN_idx; NaN];
%tr_ep = [ts(trP_idx(1:end-2))' ts(trN_idx(2:end)-1)']; % iso end-1 end-2 for previous mouse -fix this'
tr_ep = [ts(trP_idx)' ts(trN_idx)']

[status, interval] = InIntervals(ts,tr_ep); % detects what times belong to each trial

n1      = histoc(interval(interval>0),1:size(tr_ep,1)); % how many timestamps or position values should there be per bin
len_ep  = mat2cell(pos(status)',n1); % position values per trial
ts_ep   = mat2cell(ts(status)',n1); % time stamps per trial

%get smoothed position
len_ep = cellfun(@(a) nanconvn(a,k'),len_ep,'uni',0);

%get velocity
vel_ep = cellfun(@(a) diff([nan ;a]),len_ep,'uni',0);

% %find indiced where mouse runs fast, better for place tuning, otherwise
% spikes might exist not due to speed

fast_ep  =cellfun(@(a) find(-a>thres),vel_ep,'uni',0);
ts_ep_fast = cellfun(@(a,b) a(b),ts_ep,fast_ep,'uni',0);
vel_ep_fast = cellfun(@(a,b) a(b),vel_ep,fast_ep,'uni',0);
len_ep_fast = cellfun(@(a,b) a(b),len_ep,fast_ep,'uni',0);

% figure,
% plot(len_ep{1,1})
end
