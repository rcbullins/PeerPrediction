function [] = PP_Rastor_AssemblyStrength()

%4/6/20 - Reagan Bullins

lineIdx = 0;
    colorIdx = 1;
    colorsPlot = ['r','k','b','g', 'c','m'];
    figure
    count = 0;
    for ichan = 1:params.nChans
        UnitonChanIdx = find(spikes.maxWaveformCh == ichan);    
        if  isempty(UnitonChanIdx) 
            continue
        else 
            for icluster = 1: length(UnitonChanIdx) 
                lineIdx = lineIdx + 1;
                plot(spikes.times{UnitonChanIdx(icluster)}, lineIdx*ones(length(spikes.times{UnitonChanIdx(icluster)}),1), ['.' colorsPlot(colorIdx)]); 
               channel_list(lineIdx) = ichan;
               hold on
            end
            %make colors loop if go over 6
            if colorIdx < 6
                colorIdx = colorIdx + 1;
            else
                %disp("Need more colors");
                colorIdx = 1;
            end
        end
    end
    
%  
    % channel_list and color_name both give channel number and
    % corresponding color --> make legend with these somehow
    color_num = length(channel_list);
    channel_list = strsplit(num2str(channel_list));
    
    legend([channel_list], 'location', 'Best')
    
    title("Anatomical Rastor: RSC")
    xlabel("Time (s)")
    ylabel("Cell Number")
  
    xlim([xlim_lower xlim_upper]) % xline()

    hold off