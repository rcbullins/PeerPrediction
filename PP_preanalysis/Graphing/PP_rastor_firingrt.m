function PP_rastor_firingrt(params, spikes, xlim_lower, xlim_upper)

%Still needs work

     FiringRateChan = cellfun(@length, spikes.times)/params.length_recording
     %Make Matrix with firgin rate on left and cell number on right :)
     cluster_vector = (1:1:92);
     paired_FiringRate_CellNumIdx = [FiringRateChan(:), cluster_vector(:)]
     %paired_FiringRate_CellNumIdx = [FiringRateChan(:), spikes.maxWaveformCh(:)]
     %Sort Matrix to be in order of least firing to most
     paired_SortedFiringRate_CellNumIdx = sortrows(paired_FiringRate_CellNumIdx, 1);
     sorted_CellNumIdx = paired_SortedFiringRate_CellNumIdx(:,2);
    lineIdx = 0;
    colorIdx = 1;
    colorsPlot = ['k','r','b','g', 'c','m'];
    figure
    count = 0;
   
    
     for icluster = 1: length(sorted_CellNumIdx) 
         lineIdx = lineIdx + 1;
         plot(spikes.times{sorted_CellNumIdx(icluster)}, lineIdx*ones(length(spikes.times{sorted_CellNumIdx(icluster)}),1), ['.' colorsPlot(colorIdx)]); 
         colorIdx = colorIdx +1;
         hold on
         if colorIdx < 6
                colorIdx = colorIdx + 1;
         else
                colorIdx = 1;
         end
     end
          
    
    %legend([sorted_CellNumIdx], 'location', 'Best')
    
    title("Firing Rate Rastor: RSC")
    xlabel("Time (s)")
    ylabel("Cell Number")
  
    xlim([xlim_lower xlim_upper]) % xline()

    hold off
end