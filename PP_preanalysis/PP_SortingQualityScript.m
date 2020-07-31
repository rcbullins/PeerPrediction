%PP_sorting_quality_script

%% Cluster Quality
%resultsDirectory = 'C:\Users\rcbul\Documents\English Lab\PP_RSC_Data\u19_200313_155505'
resultsDirectory = 'C:\Users\rcbul\Documents\English Lab\PP_RSC_Data\m115_191203_152410_n'
[clusterIDs, unitQuality, contaminationRate] = maskedClusterQualityKilosort(resultsDirectory) %kilosort

%% Set Graph Defaults Now
SetGraphDefaults;

%% 
%hpc super pyram
super_idx = [2 4 5 6 10:20 22 23 25:32 34 35]
super_pyram = supergoodUnitQualIdx(super_idx,:)

%% Plotting All Clusters from Function
    plot(unitQuality(:,1), contaminationRate(:,1), '.b');
    title('Cluster Quality');
    xlabel('Unit Quality');
    ylabel('Contamination Rate');
    hold on
    plot_these = find(unitQuality >= 20);
    plot(unitQuality(plot_these,1), contaminationRate(plot_these,1),'.r');
    legend('Units < 20 Isolation quality','Units >= 20 Isolation quality')
    
    %Minus 1, 0 based
    clusterIDs_f = clusterIDs_f - 1;
    %%
    idx_units =  zeros(1, length(clusterIDs_p));
       
    for icluster = 1:length(clusterIDs_p) %manual sorted
       idx_units(1,icluster) = find(clusterIDs_f == clusterIDs_p(icluster)); 
    end
    
    goodUnitQual = unitQuality(idx_units(1,:));
    goodclusterIDs= clusterIDs_f(idx_units(1,:));
    supergoodUnitQualIdx = find(goodUnitQual>=20);
    
    gscount =0;
    GoodSpikes = {};
    for i = supergoodUnitQualIdx'
        gscount = gscount+1;
        GoodSpikes{gscount} = spikes.times{i};
    
    end
    plot(unitQuality(idx_units(1,:)), contaminationRate(idx_units(1,:)),'.r');
    
    