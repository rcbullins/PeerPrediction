function [] = Concat_Sessions_OptHist(winRange, num_mat_files)
% Purpose: Load and concatenate data from specified sessions, to make a
%          summary optimal time window histogram
% Input:   winRange, num_mat_files (specified num title of sessions to use)
% Output:  histogram of all data from specified sessions
% Dependencies: Data from assembly function
%               mat files created with optimal win for each session
% Created: 9/16/20 by Reagan Bullins


%% Load specified mat files 
 for file_num = 1:length(num_mat_files)
     load(['mat' num2str(num_mat_files(file_num)) '.mat'])
 end

%% Put all optimal windows in an array
 optimal_win = {}
 optimal_win{1} = optimal_win2;
 optimal_win{2} = optimal_win3;
 %optimal_win{3} = optimal_win5;
 optimal_win{3} = optimal_win6;

%% Count per session how many optimal windows per bin
    ct_per_win_mat = zeros(length(num_mat_files), length(winRange));
 
  % for every session
    for isess = 1:length(num_mat_files)
      % for every window, count number of optimal time windows
        for iwin = 1:length(winRange)
            ct_per_win_mat(isess,iwin) = sum(optimal_win{isess}(:,1) == winRange(1,iwin));
        end
    end
%% Sum columns together to get wholistic number of counts per bin
    total_ct_per_win = zeros(1, length(winRange));
    for icol = 1:length(winRange)
    total_ct_per_win(1,icol) = sum(ct_per_win_mat(:,icol));
    end
%% Graph
    bar(1:length(winRange), total_ct_per_win, 'FaceColor', 'b', 'facealpha',.5, 'edgecolor', 'none');
    hold on

    title({'Multiple Sessions';'Optimal Time Window Histogram'});
    ylabel('Count');
    xlabel('Log Scale Time Windows (ms)');
    
    %set(gca,'XTickLabel', [1 2 4 8 16 32 64 128 256 512 1024])
    set(gca,'XTickLabel', [0 16 26 36 46 56 128])
    
    optimal_win_tot = [optimal_win{1}(:,1); optimal_win{2}(:,1); optimal_win{3}(:,1)] %; optimal_win{4}(:,1)];
    med_opt1 = (['Median Window = ' num2str(median(optimal_win_tot)*1000)]);
    text(.5, max(total_ct_per_win), med_opt1, 'Color', 'b');

end