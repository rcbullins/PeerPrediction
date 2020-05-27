function [smoothedTrains, pair_idx] = PP_Raster_SmoothedTrains(spike_info_path, pairsToRun)

%% Rastor for Actual and Predictor
cd(spike_info_path)
load('m115_191203_152410_2 peerPrediction_inputs.mat')

prompt = 'Which Cells to Compare? Format in [predict_cell#, act_cell#]: '
cells_to_compare = input(prompt)

predict_cell = cells_to_compare(1);
act_cell = cells_to_compare(2);

%% Get smoothed spike trains for different smoothing windows

spikeTimes = binned_spikes(:,1,:); %FOR NOW with only ONE Predictor

for win = 1:length(winRange)
    for time_win = 1:size(spikeTimes, 3)%for all time windows
        if time_win > win & time_win + win < size(spikeTimes,3)
             if win == 0 
                smoothedTrains(:,time_win)= spikeTimes(predict_cell,:,time_win);
             elseif win > 0
                smoothedTrains(:,time_win)=  sum(spikeTimes(predict_cell,:,(time_win - win:time_win + win))/(win*2),3);
             end
        end
            if time_win <= win
                j = time_win - 1;
                smoothedTrains(:,time_win)=  sum(spikeTimes(predict_cell,:,(time_win - j:time_win + win))/(length(1:time_win + win)),3);
            end
            if time_win + win >= size(spikeTimes,3)
                j = size(spikeTimes,3) - time_win;
                smoothedTrains(:,time_win)=  sum(spikeTimes(predict_cell,:,(time_win - win:time_win + j))/(length(time_win - win:size(spikeTimes,3))),3);
            end 
    end
end
%% Find the pair index out of all pairs for the cells you want to compare
%Account for target cell may or may not be larger than the comparison cell
if predict_cell < act_cell
    row_idx = find(pairsToRun(:,1) == predict_cell)
    chosen_idx = find(pairsToRun(row_idx,2) == act_cell)
    pair_idx = row_idx(chosen_idx)
else 
    row_idx = find(pairsToRun(:,2) == predict_cell)
    chosen_idx = find(pairsToRun(row_idx,1) == act_cell)
    pair_idx = row_idx(chosen_idx)
end
%%
%Initiate Rastor Plot
%plot the target cell -- PREDICT or ACTUAL??
        subplot(6,1,1)
    for ispike = 1:length(spikes.times{predict_cell})
        xline(spikes.times{predict_cell}(ispike))
        hold on
    end
%plot compairson cell    
    subplot(6,1,2)
     for ispike = 1:length(spikes.times{act_cell})
        xline(spikes.times{act_cell}(ispike))
        hold on
     end
%plot smoothed train at 25 ms
    x = (0:150)
    subplot(6,1,3)
    for ispike = 1:length(smoothedTrains)
        plot(x, smoothedTrains(ispike, 25))
        hold on
    end

%plot smoothed train at 50 ms
    subplot(6,1,4)
    for ispike = 1:length(smoothedTrains)
        xline(smoothedTrains(ispike, 50))
        hold on
    end

%plot smoothed train at 100 ms
    subplot(6,1,5)
    for ispike = 1:length(smoothedTrains)
        xline(smoothedTrains(ispike, 100))
        hold on
    end

%plot smoothed train at 150 ms
    subplot(6,1,6)
    for ispike = 1:length(smoothedTrains)
        xline(smoothedTrains(ispike, 150))
        hold on
    end
end