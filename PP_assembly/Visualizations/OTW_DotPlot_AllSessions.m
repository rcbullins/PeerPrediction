function [] = OTW_DotPlot_AllSessions(winRange, path_mat_files1, path_mat_files2, num_mat_files)
% Purpose: Show the median optimal time window for each session as a
%           different colored dot. Also show the median optimal window for
%           all sessions as a line. Additionally compare all cells vs only
%           quality cells optimal median time windows.
% Input:   winRange, num_mat_files (specified num title of sessions to use)
%           path_mat_files1 (all cells)
%           path_mat_files2 (quality cells)
% Output:  Dot plot, each dot represents a different session
% Dependencies: Data from assembly function
%               mat files created with optimal win for each session
% Created: 10/23/20 by Reagan Bullins
 
%% First load in data for first session set
cd(path_mat_files1)
    %% Load specified mat files 
     for file_num = 1:length(num_mat_files)
         load(['mat' num2str(num_mat_files(file_num)) '.mat'])
     end

    %% Put all optimal windows in an array
     optimal_win_1 = zeros(1,length(num_mat_files))
     optimal_win_1(1,1) = median(optimal_win1);
     optimal_win_1(1,2) = median(optimal_win2);
     optimal_win_1(1,3) = median(optimal_win3);
     optimal_win_1(1,4) = median(optimal_win5);
     optimal_win_1(1,5) = median(optimal_win6);
%% Load in data for second session set
clear optimal_win1 optimal_win2 optmial_win3 optimal_win5 optimal_win6
cd(path_mat_files2)
    %% Load specified mat files 
     for file_num = 1:length(num_mat_files)
         load(['mat' num2str(num_mat_files(file_num)) '.mat'])
     end

    %% Put all optimal windows in an array
     optimal_win_2 = zeros(1,length(num_mat_files))
     optimal_win_2(1,1) = median(optimal_win1);
     optimal_win_2(1,2) = median(optimal_win2);
     optimal_win_2(1,3) = median(optimal_win3);
     optimal_win_2(1,4) = median(optimal_win5);
     optimal_win_2(1,5) = median(optimal_win6);
     clear optimal_win1 optimal_win2 optimal_win3 optimal_win5 optimal_win6
%% Plot Data
    color_dot = ['.r';'.b';'.g';'.y';'.c';'.p']
    figure
    %small bug -- if median ends up being odd this won't work
    for i = 1:length(num_mat_files)
        if isnan(optimal_win_1(1,i)) == 0   
            disp('why')
            winRange_idx = find(optimal_win_1(1,i) == winRange(1,:))
            plot(1, winRange_idx, color_dot(i,:), 'MarkerSize',15)
            hold on
        end
        
        if isnan(optimal_win_2(1,i)) == 0
            winRange_idx2 = find(optimal_win_2(1,i) == winRange(1,:))
            plot(1.5, winRange_idx2, color_dot(i,:), 'MarkerSize',15)
            hold on
        end
    end
    %this is all hard coded
     set(gca,'YTickLabel', [winRange(18)*1000 winRange(20)*1000 winRange(22)*1000 ...
         winRange(24)*1000 winRange(26)*1000 winRange(28)*1000 winRange(30)*1000 winRange(32)*1000 1.024*1000])
    xlim([.5 2])
    set(gca, 'XTickLabel',[1 1.5])
    xticklabels({'','All Cells','Quality Cells'})
    title({'Optimal Median Time Window';'per session'})
    ylabel('Median Optimal Time Window')
    
    %% add in median lines
    %one option median of median
    optimal_win_1(find(isnan(optimal_win_1))) = [];
    median1 = find(median(optimal_win_1(1,:)) == winRange(1,:))
    plot([.9 1.1], [median1 median1], '-k')
  
    optimal_win_2(find(isnan(optimal_win_2))) = []
    median2 = find(median(optimal_win_2(1,:)) == winRange(1,:))
    plot([1.4 1.6], [median2 median2], '-k')
    

end