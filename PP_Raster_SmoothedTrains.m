function [smoothedTrains, pair_idx] = PP_Raster_SmoothedTrains(choice_graphPairs, binned_spikes, pairsToRun, spikes)
%Purpose: Create cross correlograms for specified pairs.

%Dependencies: Buzcode
%              Paths set in beggining of PP_Analyzing_Script
%              Spiking Data from Recording

%Inputs: choice_graphPairs (what pairs to graph)
%        binned_spikes (from pp_inputs)
%        

%Output: Creates a figure with 6 subplots for each pair given
%           First Raster: Cell 1 Spikes (actual)
%           Second Raster: Cell 2 Spikes (predictor)
%           Third-sixth Raster: Cell 2 Spikes smoothed over varying windows
%        smoothedTrains: smoothed spike times of cell 2 (only specified
%           graphing windows, or you can change code to run all)
%        pairIdx: index of pair in dev

%Created: 06/03/20 by Reagan Bullins

%% Get smoothed spike trains for different smoothing windows

spikeTimes = binned_spikes(:,1,:); %FOR NOW with only ONE Predictor
%winRange = (0:150); %if want all windows smoothed over, & change
%   win_num(win) to win 
win_num = [25 50] % 100 150]; %the ms windows we will look at FOR SPEED
for ipair = 1:length(choice_graphPairs)
    figure
    %identify which number cells are being compaired in this pair
    cell_1 = pairsToRun(choice_graphPairs(ipair), 1);
    cell_2 = pairsToRun(choice_graphPairs(ipair), 2);
    %Get the spike times for these two cells
    actual_cell = spikes.times{cell_1};
    predict_cell = spikes.times{cell_2};
  %Smooth predict cell
        for win = 1:length(win_num)
           disp( 'hi');
            for time_win = 1:size(spikeTimes, 3)%for length of whole recording
                if time_win > win_num(win) & time_win + win_num(win) < size(spikeTimes,3)
                     if win_num(win) == 0 
                        smoothedTrains(:,time_win)= spikeTimes(cell_2,:,time_win);
                     elseif win_num(win) > 0
                        smoothedTrains(:,time_win)=  sum(spikeTimes(cell_2,:,(time_win - win_num(win):time_win + win_num(win)))/(win_num(win)*2),3);
                     end
                end
                    if time_win <= win_num(win)
                        j = time_win - 1;
                        smoothedTrains(:,time_win)=  sum(spikeTimes(cell_2,:,(time_win - j:time_win + win_num(win)))/(length(1:time_win + win_num(win))),3);
                    end
                    if time_win + win_num(win) >= size(spikeTimes,3)
                        j = size(spikeTimes,3) - time_win;
                        smoothedTrains(:,time_win)=  sum(spikeTimes(cell_2,:,(time_win - win_num(win):time_win + j))/(length(time_win - win_num(win):size(spikeTimes,3))),3);
                    end 
            end
        end
%================================HERE
    %Initiate Rastor Plot
%plot the target cell -- ACTUAL cell
        subplot(6,1,1)
        plot(actual_cell, ones(length(actual_cell)),'.r')
        subplot(6,1,2)
        plot(predict_cell, ones(length(predict_cell)),'.r')
        
        subplot(6,1,3)
       % plot(smoothedTrains(:, 25),ones(length(smoothedTrains(:)),'.r')
      
        subplot(6,1,4)
        plot(smoothedTrains(:, 50),'.r')
       
        subplot(6,1,5)
        plot(smoothedTrains(:, 100)'.r')
        
        subplot(6,1,6)
        plot(smoothedTrains(:, 150),'.r')
        
% ===================================Old        
% %plot compairson cell    
%     subplot(6,1,2)
%      for ispike = 1:length(predict_cell)
%         plot(predict_cell(ispike), '.')
%         hold on
%      end
% %plot smoothed train at 25 ms
%     x = (0:150)
%     
%     for ispike = 1:length(smoothedTrains)
%         if ispike == 10
%             disp('well hi')
%         end
%         subplot(6,1,3)
%         plot(x, smoothedTrains(ispike, 25),'.')
%         hold on
%         subplot(6,1,4)
%         plot(smoothedTrains(ispike, 50),'.')
%         hold on
%         subplot(6,1,5)
%         plot(smoothedTrains(ispike, 100)'.')
%         hold on
%         subplot(6,1,6)
%         plot(smoothedTrains(ispike, 150),'.')
%         hold on
%     end

    disp('Yay Science!')
end