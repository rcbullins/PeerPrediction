function [smoothedTrains] = PP_Raster_SmoothedTrains(choice_graphPairs,choice_sec, binned_spikes, pairsToRun, spikes)
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
%Options: win_num: windows to smooth data over, you can change this

%Credits: Adaptations to David Tingley's Smooth Train Code is implemented
%Created: 06/03/20 by Reagan Bullins

%% Get smoothed spike trains for different smoothing windows

spikeTimes = binned_spikes(:,1,:); %FOR NOW with only ONE Predictor
%winRange = (0:150); %if want all windows smoothed over, & change
%   win_num(win) to win 
win_num = [25 50 100 150]; %the ms windows we will look at FOR SPEED

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
                        smoothedTrains(win,time_win)= spikeTimes(cell_2,:,time_win);
                     elseif win_num(win) > 0
                        smoothedTrains(win,time_win)=  sum(spikeTimes(cell_2,:,(time_win - win_num(win):time_win + win_num(win)))/(win_num(win)*2),3);
                     end
                end
                    if time_win <= win_num(win)
                      
                        j = time_win - 1;
                        smoothedTrains(win,time_win)=  sum(spikeTimes(cell_2,:,(time_win - j:time_win + win_num(win)))/(length(1:time_win + win_num(win))),3);
                    end
                    if time_win + win_num(win) >= size(spikeTimes,3)
                      
                        j = size(spikeTimes,3) - time_win;
                        smoothedTrains(win,time_win)=  sum(spikeTimes(cell_2,:,(time_win - win_num(win):time_win + j))/(length(time_win - win_num(win):size(spikeTimes,3))),3);
                    end 
            end
        end
%================================
%Graphing
    %plot ACTUAL cell
        subplot(6,1,1);
        title('Actual Cell');
        %find the index points of the second you want to graph
        actual_x = find(actual_cell >= choice_sec & actual_cell < choice_sec+1);
        for idx_x = 1:length(actual_x)
            xline(actual_cell(actual_x(idx_x)));
        end
         %plot(actual_cell(actual_x),ones(length(actual_x),'.r'));
         xlim([choice_sec choice_sec+1]);
         set(gca,'XTick',[]);
         set(gca,'YTick',[]);
    %plot PREDICT cell   
        subplot(6,1,2)
        title('Predictor Cell')
        predict_x = find(predict_cell >= choice_sec & predict_cell < choice_sec+1);
        for idx_x = 1:length(predict_x)
            xline(predict_cell(predict_x(idx_x)));
        end
        %plot(predict_cell(predict_x),ones(length(predict_x),'.r'));
        xlim([choice_sec choice_sec+1]);
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);
    %plot Smoothing windows
        link_ax = zeros(length(win_num));
        for idx_win = 1:length(win_num)
            link_ax(idx_win) = subplot(6,1,idx_win+2);
            plot(smoothedTrains(idx_win, choice_sec*1000:(choice_sec+1)*1000));
            hold on
             %sgtitle(['Actual vs Predictor Spike Trains: Cell ' num2str(cell_1) ' vs ' num2str(cell_2)']);
            title(['Smoothed Predictor ' num2str(win_num(idx_win)) 'ms'])
            xlim([0 1000])
            if idx_win < length(win_num)
                set(gca,'XTick',[]);
            elseif idx_win == length(win_num)
                 xticks([0 500 1000]);
                 xticklabels({[choice_sec, choice_sec+.5, choice_sec+1]});
                  xlabel('Time (s)');
                   ylabel('Probability');
            end
          
         end
    %Graph specs
       % sgtitle(['Actual vs Predictor Spike Trains: Cell ' num2str(cell_1) ' vs ' num2str(cell_2)']);
     
        hold off
end
    disp('Yay Science!')
end