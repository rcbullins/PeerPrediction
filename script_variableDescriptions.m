%What variables are & some code
   load('epoch_variables.mat')
%% Heat Maps: All Cells
% spk_center_epoch  : matrix (cell X epoch) with # of spikes, for center 100 ms of epoch
% spk_whole_epoch   : matrix (cell X epoch) with # spikes, for whole epoch
    figure
    imagesc(spk_center_epoch)
    title('Center Epoch: All cells')
    xlabel('Epoch')
    ylabel('Cell #')
    colorbar
  
    figure 
    imagesc(spk_whole_epoch)
    title('Whole Epoch: All cells')
    xlabel('Epoch')
    ylabel('Cell #')
    colorbar;
   
%% Heat Maps: Cells with good isolation quality
% supergoodUnitQualIdx : cells indicies with isolation quality >= 20
% spk_whole_quality: (cell X epoch) for cells with isolation quality >=20, the number of spike times per epoch
% spk_center_quality: (cell X epoch) for cells with isolation quality >= 20
    figure
    imagesc(spk_whole_quality)
    title('Whole Epoch: Good Isolation Quality')
    xlabel('Epoch')
    ylabel('Cell #')
    lcolorbar('Number of Spikes')
    figure 
    imagesc(spk_center_quality)
    title('Center Epoch: Good Isolation Quality')
    xlabel('Epoch')
    ylabel('Cell #')
    lcolorbar('Number of Spikes')
    
%% Heat Maps: Good Isolation Quality & >100 epochs
% supergoodUnitQualIdx : cells indicies with isolation quality >= 20
% whole_great_cells : cells with good iso quality and fire in >100 epochs
% center_great_cells : cells with good iso quality and fire in >100 epochs

% whole_spk_great : (cell X epoch) number of spikes for each cell with good
% iso quality and spiking in over 100 epochs
% center_spk_great: (cell X epoch) same as above
    figure
    imagesc(whole_spk_great)
    title('Whole Epoch: Good Iso Quality & > 100 Epoch')
    xlabel('Epoch')
    ylabel('Cell #')
    lcolorbar('Number of Spikes')
    colorbar
    figure 
    imagesc(center_spk_great)
    title('Center Epoch: Good Iso Quality & >100 Epoch')
    xlabel('Epoch')
    ylabel('Cell #')
    lcolorbar('Number of Spikes')
    colorbar

%% How many epochs each cell participates in

% whole_epoch_sum  : (cell X 1) for each cell, how many epochs does it fire in
% center_epoch_sum : (cell X 1) for each cell, how many center epochs does it fire in

% whole_quality_sum : (cell X 1) for each quality cell (>=20), how many epochs does it fire in
% center_quality_sum: (cell X 1) for each quality cell (>=20), how many center epochs does it fire in

%% new spikes to use

spikeGreatTimes = {};

for icell = 1:length(center_great_cells)
   spikeGreatTimes{icell}(:) = spikes.times{center_great_cells(icell)}(:);
end

epochPulses = pulseEpochs;
epochPulses(:,1) = pulseEpochs(:,1) +.1;
epochPulses(:,2) = pulseEpochs(:,2) - .1;




