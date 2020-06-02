function [smoothedTrains, pair_idx] = PP_Raster_SmoothedTrains(choice_graphPairs, spike_info_path, basepath)

%% Rastor for Actual and Predictor
cd(spike_info_path)
load('m115_191203_152410_2 peerPrediction_inputs.mat')
load('m115_191203_152410_2.spikes.cellinfo.mat')

cd(basepath)
load('pairsToRun.mat')

%% Get smoothed spike trains for different smoothing windows

spikeTimes = binned_spikes(:,1,:); %FOR NOW with only ONE Predictor

%% Find the pair index out of all pairs for the cells you want to compar

for ipair = 1:length(choice_graphPairs)
    figure
    %identify which number cells are being compaired in this pair
    cell_1 = pairsToRun(choice_graphPairs(ipair), 1);
    cell_2 = pairsToRun(choice_graphPairs(ipair), 2);
    %Get the spike times for these two cells
    actual_cell = spikes.times{cell_1};
    predict_cell = spikes.times{cell_2};
  %Smooth predict cell
        for win = 1:length(winRange)
            for time_win = 1:size(spikeTimes, 3)%for all time windows
                if time_win > win & time_win + win < size(spikeTimes,3)
                     if win == 0 
                        smoothedTrains(:,time_win)= spikeTimes(cell_2,:,time_win);
                     elseif win > 0
                        smoothedTrains(:,time_win)=  sum(spikeTimes(cell_2,:,(time_win - win:time_win + win))/(win*2),3);
                     end
                end
                    if time_win <= win
                        j = time_win - 1;
                        smoothedTrains(:,time_win)=  sum(spikeTimes(cell_2,:,(time_win - j:time_win + win))/(length(1:time_win + win)),3);
                    end
                    if time_win + win >= size(spikeTimes,3)
                        j = size(spikeTimes,3) - time_win;
                        smoothedTrains(:,time_win)=  sum(spikeTimes(cell_2,:,(time_win - win:time_win + j))/(length(time_win - win:size(spikeTimes,3))),3);
                    end 
            end
        end

    %Initiate Rastor Plot
%plot the target cell -- ACTUAL cell
        subplot(6,1,1)
    for ispike = 1:length(actual_cell)
        xline(actual_cell(ispike))
        hold on
    end
%plot compairson cell    
    subplot(6,1,2)
     for ispike = 1:length(predict_cell)
        xline(predict_cell(ispike))
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
    disp('Yay Science!')
end