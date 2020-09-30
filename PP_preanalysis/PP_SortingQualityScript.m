%PP_sorting_quality_script
%Quality and Inclusion Criteria -(firing rate)

%% Cluster Quality
%resultsDirectory = 'C:\Users\rcbul\Documents\English Lab\PP_RSC_Data\u19_200313_155505'
resultsDirectory = 'C:\Users\rcbul\Documents\English Lab\PP_RSC_Data\m115_191203_152410_n'
[clusterIDs, unitQuality, contaminationRate] = maskedClusterQualityKilosort(resultsDirectory) %kilosort

%% Set Graph Defaults Now
SetGraphDefaults;

%% Plotting All Clusters from Function
    plot(unitQuality(:,1), contaminationRate(:,1), '.b');
    title('Cluster Quality');
    xlabel('Unit Quality');
    ylabel('Contamination Rate');
    hold on
    plot_these = find(unitQuality >= 20);
    plot(unitQuality(plot_these,1), contaminationRate(plot_these,1),'.r');
    legend('Units < 20 Isolation quality','Units >= 20 Isolation quality')
    
   
    %% Find out what clusters are quality
     %Minus 1, 0 based
    clusterIDs_f = clusterIDs - 1;
  
    clusterIDs_p = spikes.UID;
    idx_units =  zeros(1, length(clusterIDs_p));
    
    %find clusters from kilosort and manual sorted, if do not exist make it
    %a NaN value
    for icluster = 1:length(clusterIDs_p) %manual sorted
        exist_cluster = find(clusterIDs_f == clusterIDs_p(icluster))
        if ~isempty(exist_cluster)
            idx_units(1,icluster) = find(clusterIDs_f == clusterIDs_p(icluster)); 
        else
             idx_units(1,icluster) = NaN; 
        end
    end
    %Take out NaN values
    idx_units = idx_units(~isnan(idx_units(1,:)));
    %Find good quality clusters, then find the index for ones with quality
    %over 20
    goodUnitQual = unitQuality(idx_units(1,:));
    goodclusterIDs= clusterIDs_f(idx_units(1,:)); %units on both manual and computer sorted
    supergoodUnitQualIdx = find(goodUnitQual>=20); %units quality greater than 20 - is the idx in spikes struct

 %% Inclusion Criteria Pt II - Firing Rate of Cells
 load([session_name '.spikes.cellinfo.mat']);
 firing_rate = zeros(length(spikes.times),1);
    for icell = 1:length(spikes.times)
    num_spk = length(spikes.times{icell});
    length_time = spikes.times{icell}(length(spikes.times{icell})) - spikes.times{icell}(1);
    firing_rate(icell,1) = num_spk/length_time;
    end
 % need a threshold for what firing rate is a good one
 %QUESTIION -- this does not mean consistency, how can i figure this out
cell_idx_good_fr = find(firing_rate(:,1) > 1) %LOOK INTO THIS 
    
%% Now see what cells have BOTH good quality and consistent firing rate