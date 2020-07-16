function [] = Weighted_Raster_FiringRate(spikes, weights, target_cell, time_plot, optimal_win, winRange)
%Purpose: Creates a raster plot for a target cell. The target cell's raster
    %will be at the top and the second half of the figure will contain peer
    %cell raster plots ordered on the y-axis by weights. Below the black line
    %will be negative weighted cell activity, and above the black line will be
    %positive weighted cell activity. The color bar depicts the firing rate
    %of each peer cell. 
%Inputs: spikes (spikes struct)
%        weights (given by CrossValidationAssemblyPrediction function)
%        target cell (defined by user)
%        time_polt (this is the actual second of the recording to be plotted)
%        optimal_win (optimal window per cell)
%        winRange (bin sizes)
%Outputs: raster plot ordered by magnituded of weight of peer cell
%Dependencies: crossValidationAssemblyPrediction

%Created 7/16/20 by Reagan Bullins (variation of weighted_raster_Asmb)

%%
% Find optimal window for target_cell
    tar_opt_win = optimal_win(target_cell,1);
% Find which index in winRange this window is
    win_idx = find(winRange == tar_opt_win)
% Get Weights of peer cells for target cell
    target_weights = (weights{win_idx}(target_cell+1,:)); 
    num_cells = size(target_weights, 2)
    idx_weights = (1:num_cells);
% Take out nan value (given target cell idx)
    idx_nan = find(isnan(target_weights));
    target_weights(idx_nan) = [];
    idx_weights(idx_nan) = [];
    zero_index = find(target_weights < 0, 1, 'first');
% Create matrix of indexes and corresponding weights for each peer cell
    target_idx_mat = [idx_weights;target_weights]';
% Sort the peer cells based on weights (most - to most +)
    sorted_weights = sortrows(target_idx_mat,2);
% Add in row of zeros where weights change negative to positive: graphing
% purposes
    zero_index = find(sorted_weights(:,2) > 0, 1, 'first');
    sorted_weights_wZ = zeros(num_cells,2);
    sorted_weights_wZ(:,1) = [sorted_weights(1:zero_index-1,1)' 0 sorted_weights(zero_index:length(sorted_weights),1)']';
    sorted_weights_wZ(:,2) = [sorted_weights(1:zero_index-1,2)' 0 sorted_weights(zero_index:length(sorted_weights),2)']';
% Get index for raster plotting
    raster_idx = sorted_weights_wZ(:,1);

%% plot target cell
    h1 = subplot(2,1,1);
    title(['Target Cell: ' num2str(target_cell)])
    target_spikes = spikes.times{target_cell};
   % target_x = find(target_spikes >= time_plot & target_spikes < time_plot+1);
        for idx_x = 1:length(target_spikes)
            %xline(target_spikes(target_x(idx_x)), 'LineWidth', 1.25, 'Color','k');
            xline(target_spikes(idx_x), 'LineWidth', 1.25, 'Color','k');
        end
         xlim([time_plot time_plot+1]);
         set(gca,'XTick',[]);
         set(gca,'YTick',[]);
         
%% Get firing rate

    %want to order colors by firing rate
    firing_rate = zeros(length(sorted_weights), 2);
    max_time = max(cellfun(@max,spikes.times));
    for icell = 1:length(sorted_weights)
        num_spikes = length(spikes.times{sorted_weights(icell,1)});
        firing_rate(icell, 2) = num_spikes/max_time;
        firing_rate(icell, 1) = sorted_weights(icell, 1);
    end
    sorted_FR = sortrows(firing_rate,2);
    %cell number ordered by firing rate
    idx_FR = sorted_FR(:,1);
    colors = jet(num_cells);
    
    
%% make raster plot
    h2 = subplot(2,1,2)
  
    % for each cell, plot a raster
    for idx_cell = 1:length(raster_idx)
       % if raster_idx is ZERO, this is a marker on the graph
       if raster_idx(idx_cell) == 0
           yline(idx_cell, '--');
       else
       %current cell spikes (positive to negative weight)
       peer_cell = spikes.times{raster_idx(idx_cell)};
       % only plot given window to plot
         %predict_x = find(peer_cell >= time_plot & peer_cell < time_plot+1);
       % get y-axis established
          y_idx = (1:length(peer_cell)); %changed to peer cell from predict_x
       % which y value to plot this cell on
       y_idx(:) = idx_cell;
       %plot(peer_cell(predict_x),y_idx, '.r');
       
       %what index of firing rate this cell is, and assign color based on
       %that
       color_idx = find(idx_FR == raster_idx(idx_cell));
       scatter(peer_cell, y_idx, 4, colors(color_idx,:), 'filled');
       %xline(peer_cell(predict_x))
       hold on
       xlim([time_plot time_plot+1]);
    end
    end 
    %% Finishing touches
     xlabel('Time(s)')
     ylabel('Peer Cells Ordered by Weights')
     set(gca, 'YTick',[])
     %set positions of subplots
     set(h1, 'OuterPosition',[0,0.85,.87,.1]);
     set(h2, 'OuterPosition',[0,.1,1,.75]);
     % add color bar without y ticks
     c = colorbar('YTick', []);
     c.Ticks = [.5];
     c.TickLabels = {'Firing Rate'};
     c.Color = [0 0 0];
     text(time_plot+1.15, length(raster_idx), 'High', 'Color', colors(length(colors),:));
     text(time_plot+1.15, 1, 'Low','Color', colors(1,:));
     %make colorbar same colors as on plot
     colormap(jet);
     end