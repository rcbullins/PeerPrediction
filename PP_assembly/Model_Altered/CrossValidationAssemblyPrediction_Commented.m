function [L,weights, L_V] = CrossValidationAssemblyPrediction_ExtraPredict(spikes,velocities, varargin)
warning off
%

%  spikes = bz_code standard
%  spikes = bz_GetSpikes('basepath',basepath);
%  [L,weights] = CrossValidationAssemblyPrediction_Cells(spikes)
% https://github.com/buzsakilab/buzcode
%
% computes the loglikelihood prediction of the spike train Starget from the
% binned spike train matrix of its peers (Qassembly). The prediction is built
% upon a Generalized Linear Model of the multivariate data given by
% Qassembly. The algorithm searches for a set of weights that maximize the
% prediction of the binned spike train Qtarget (the training set) and then
% computes the likelihood (during the test set). Everything is the same as
% in Harris et al., Nature, 2004.
%
%
% Dependencies: Epoch2Mat
%               ComputePeerPredictionRoot
%               SpkTrainLogLikelihood
%               modifiedExp


% Update by RB on 7/10/20
%   -adding in velocity as an extrapredictor
%
%% default parameters
p = inputParser;
addParameter(p,'dt',.2,@isnumeric)
addParameter(p,'std',10,@isnumeric)
addParameter(p,'epoch',[0 inf],@isnumeric)
addParameter(p,'cells',[],@isnumeric)
addParameter(p,'nbEp',10,@isnumeric)
parse(p,varargin{:})
dt = p.Results.dt; %time bins
std = p.Results.std;
epoch = p.Results.epoch;
cells = p.Results.cells;
nbEp = p.Results.nbEp;

%hard coded parameters for testing -RB
% dt = .2
% std = 10
% epoch = [0 inf]
% %cells = []
% nbEp = 10 %test sets

%% Get Velocity -- working on 7/8/20 RB
    ts_ex = 0:dt*1000:length(velocities)+dt*1000; 
    [ts1_ex] = Epoch2Mat(ts_ex,epoch); %length of epoch in dt binning increments
     ts1_ex = cell2mat(ts1_ex');
     crossBins_ex = floor(linspace(1,length(ts1_ex),nbEp+1));
     crix_ex = [crossBins_ex(1) crossBins_ex(2);crossBins_ex(2:end-1)'+1 crossBins_ex(3:end)'];
     
    %average velocity data in 200 ms bins (or dt size)
    vel_avg = zeros(length(ts1_ex)-1, 1);
    %for every dt interval get an average velocity
    for i_interval = 1:length(ts1_ex)-1
        if i_interval ~= length(ts1_ex)-1
            vel_avg(i_interval,:) = mean(velocities(1, ts1_ex(i_interval)+1:ts1_ex(i_interval+1)));
        else 
            vel_avg(i_interval,:) = mean(velocities(1, ts1_ex(i_interval)+1:size(velocities, 2)));
        end
    end
   
    
    
%%
k = gaussian2Dfilter([std*10 1],[std 1]);

fprintf('Launching Cross-validated Peer info\n')


nbC = length(spikes.times); %number of cells

L = NaN(nbEp,nbC); %number of epochs x number of cells


weights = nan(nbC+1,nbEp,nbC); %initializing matrix of weights

maxT = max(cellfun(@max,spikes.times))+dt; %finding max time


ts = 0:dt:maxT; %binning time

[ts1] = Epoch2Mat(ts,epoch); %length of epoch in dt binning increments
ts1 = cell2mat(ts1');
crossBins = floor(linspace(1,length(ts1),nbEp+1));
crix = [crossBins(1) crossBins(2);crossBins(2:end-1)'+1 crossBins(3:end)'];
% for every epoch get weights
for ii=1:nbEp
    fprintf('.')
    
  %//////////RB ADD 7/10/20
     %Define velocity for this epoch
     if ii ~= nbEp
        vel_ep = vel_avg(crix(ii,1):crix(ii,2));
     else %if last one, due to averaging of velocity, have one less value
         vel_ep = vel_avg(crix(ii,1):crix(ii,2)-1);
         %add in a zero to match matrix deminsions of weights :) ?? okay??
         vel_ep(length(vel_ep)+1,1) = 0
     end
  %//////////////
  
    ep = crix(~ismember(crix,crix(ii,:),'rows'),:);
    
    bints = cellfun(@(a) ts1(a(1):a(2)),num2cell(ep,2),'uni',0); %bin edges of timestamps
    N = sum(cellfun(@numel,bints)); %adding all bin edges together
    pop = nan(N,nbC); %bins x number of cells matrix of NaN
    
    % qpop = nan(N,nbC);
    %get population
    
    %for every cell, find the timestamps in the ii epoch
    for jj = 1:nbC
        spk = cell2mat(Epoch2Mat(spikes.times{jj},ts1(ep))); %spike times in first epoch of first cell
        
        pop(:,jj) = cell2mat(cellfun(@(a) histoc(spk,a),bints,'uni',0)); %number of spikes per bin
        % pop(:,jj) = cell2mat(cellfun(@(a) nanconvn(histoc(spk,a),k),bints,'uni',0));
    end
    
    
    % for every cell, figure out what bins have spikes in them
    for kk = 1:nbC
        q = pop(:,kk);
        oc = ~ismember(1:nbC,kk);
        
        %if any of the bins has a spike in it (num spikes > 0)
        if any(q>0)
            
            warning off
            %if this is the first epoch
            if ii==1
                temp = nan(nbC+1,1);
                %determine weights
                w = ComputePeerPredictionRoot_Commented(q,[ones(size(pop,1),1) pop(:,oc)],dt); %weights
                
                temp([true oc]) = w;
                weights(:,ii,kk) = temp;
                % any epochs past the first
            else
                %Use previous weights as initial conditions to improve speed
                w0 = nanmean(weights([true oc],1:ii-1,kk),2);
                if ~any(isnan(w0))
                    w = ComputePeerPredictionRoot_Commented(q,[ones(size(pop,1),1) pop(:,oc)],dt,w0);
                else
                    w = ComputePeerPredictionRoot_Commented(q,[ones(size(pop,1),1) pop(:,oc)],dt);
                end
                
                temp = nan(nbC+1,1);
                
                temp([true oc]) = w;
                weights(:,ii,kk) = temp;
            end
            warning on
            
            
            %Now compute prediction
            ep  = crix(ii,:);
            
            bints = ts1(ep(1):ep(2)); %bin ends for timestamps for one epoch
            
            N = length(bints);
            
            popt = nan(N,nbC); % Initiate matrix of number of bins x number of cells
            %  qpopt = nan(N,nbC);
            
            %for every cell
            for jj = 1:nbC
                spk = cell2mat(Epoch2Mat(spikes.times{jj},ts1(ep)')); %all timestamps of one cell in a single epoch
                %  popt(:,jj) = nanconvn(histc(spk,bints),k);
                popt(:,jj) = histc(spk,bints); %number of spikes per bin for jjth cell
            end
            
            
            qt = popt(:,kk); %the number of spikes per bin for the kkth cell
            r = sum(qt)/length(qt)/dt; %average spikes per bin for the kkth cell
       %%/////////////////////////////////////////////////////////////////////////////
            %find the likelihood of spikes qt, with modified
            %exponential intensity function
            %PEER ACTIVITY alone, EQUATION 12 in Supplementary information
            LogFun = @(x) SpkTrainLogLikelihood(qt, modifiedExp([ones(size(popt,1),1) popt(:,oc)]*x),dt) - 0.25*x'*x;
            Lf = LogFun(w);
            
            %For Velocity and Peer -RB
            
            LogFun_V = @(x) SpkTrainLogLikelihood(qt, exp(vel_ep(:,1)) .* modifiedExp([ones(size(popt,1),1) popt(:,oc)]*x),dt) - 0.25*x'*x;
            Lf_V = LogFun_V(w);
            
            %For Space and Peer Activity (peer function x phase/place)
            %PLACEHOLDER function, something like this
            %LogFun_SpacePeer = @(x) SpkTrainLogLikelihood(qt, exp(x(1))*exp(x(2))+modifiedExp([ones(size(popt,1),1) popt(:,oc)]*x(3:end)),dt) - 0.25*x'*x;
            
       %//////////////////////////////////////////////////////////////////////////////////
            
            %find likelihood for average number of spikes for all cells (control)
            L0 = SpkTrainLogLikelihood(qt,repmat(r,[size(popt(:,oc),1) 1]),dt);
            
            %             disp(Lf-L0)
            %Get log likelihood by subtracting control likelihood from
            %a given cell log likelihood: do this for every iith epoch and
            %kkth cell
            L(ii,kk) = Lf-L0;
            L_V(ii,kk) = Lf_V-L0; %RB
        end
    end
end

fprintf('done\n')
totlen = length(ts1)*dt; %length of epoch in dt increments X dt (bin step)


L = nansum(L); %get rid of all extra NaN values
L = L/totlen; % Log Likelikhood for each cell / (length of the epoch in dt increments X dt)
L = L(:);

%/////RB
L_V = nansum(L_V); %get rid of all extra NaN values
L_V = L_V/totlen; % Log Likelikhood for each cell / (length of the epoch in dt increments X dt)
L_V = L_V(:);
%////////

weights = squeeze(nanmean(weights,2)); %take away all NaN values on edges & merge all epochs

%L_over_time = SpkTrainLogLikelihood_alt(q,modifiedExp(pop*x),dt);


end