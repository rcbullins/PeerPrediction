function [] = PP_OptimalWindow_ScatterPlot(min_win,dev,devControl, y_plot_limit)
%Purpose: Plot one point for each cell pair. This point will represent the
%         optimal time window (minumum deviance) for each cell pair. The y-axis will
%         be the error for each optimal time window (devControl_avg - dev).

%Dependencies: PP_DevianceGraphs

%Inputs: min_win (optimal window index, where the minimum deviance is)
%        dev & devControl (bz_PeerPrediction)
%        y_plot_limit (set y axis max bound)

%Output: ScatterPlot of Optimal Time windows vs Error for every Pair

%Created: 4/9/20 by Reagan Bullins

%%
%[bin_win_max,num_pairs] = size(dev); % CAN adjust
[bin_win_max] = size(dev, 1);
num_pairs = length(min_win);
devControl_avg = (mean(devControl(:,:,:),3));

for ipair = 1:num_pairs
    idx_dev = min_win(ipair); %min_win gives the index of where the minimum deviance is
    error = abs(devControl_avg(idx_dev, ipair) - dev(idx_dev, ipair));
    plot(min_win(ipair), error, '.')
    hold on
end
xlim([0 bin_win_max])
ylim([0 y_plot_limit])
title("Optimal Temporal Window: All Pairs")
xlabel("Smoothing Window (ms)")
ylabel("Error Difference")

end