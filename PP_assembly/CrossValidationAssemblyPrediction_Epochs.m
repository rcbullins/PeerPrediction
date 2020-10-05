function [L,weights] = CrossValidationAssemblyPrediction_Epochs(spikes,varargin)
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


%%
%k = gaussian2Dfilter([std*10 1],[std 1]);
k = imgaussfilt([std*10 1],[std 1]);
fprintf('Launching Cross-validated Peer info\n')


nbC = length(spikes.times); %number of cells

L = NaN(nbEp,nbC); %number cross validation folds by number of cells


weights = nan(nbC+1,nbEp,nbC); %initializing matrix of weights

maxT = max(cellfun(@max,spikes.times))+dt; %finding max time

ts = 0:dt:maxT; %binning time 

[ts1] = Epoch2Mat(ts,epoch); %what ts are in the epoch given
%-----------Reagan Editing Here ------------%
%Taking out ' from original fixes it....?
    %ts1 = cell2mat(ts1'); %original Sams
    ts1 = cell2mat(ts1); %gives consecutive timestamps for all epochs
%-------------------------------------------%

crossBins = floor(linspace(1,length(ts1),nbEp+1)); %
crix = [crossBins(1) crossBins(2);crossBins(2:end-1)'+1 crossBins(3:end)'];
for ii=1:nbEp %how many times it cross validates? nbEp = 10 folds?
    fprintf('.')
    
    ep = crix(~ismember(crix,crix(ii,:),'rows'),:);
    
    bints = cellfun(@(a) ts1(a(1):a(2)),num2cell(ep,2),'uni',0);
    N = sum(cellfun(@numel,bints));
    pop = nan(N,nbC); %bin matrix
   % qpop = nan(N,nbC);
    %get population
    
    for jj = 1:nbC
        spk = cell2mat(Epoch2Mat(spikes.times{jj},ts1(ep)));
        
        pop(:,jj) = cell2mat(cellfun(@(a) histoc(spk,a),bints,'uni',0));
       % pop(:,jj) = cell2mat(cellfun(@(a) nanconvn(histoc(spk,a),k),bints,'uni',0));
        
    end
    
    for kk = 1:nbC
        q = pop(:,kk);
        oc = ~ismember(1:nbC,kk);
        if any(q>0)
            
            warning off
            if ii==1
                temp = nan(nbC+1,1);
                w = ComputePeerPredictionRoot(q,[ones(size(pop,1),1) pop(:,oc)],dt); %weights
               
                temp([true oc]) = w;
                weights(:,ii,kk) = temp;
            else
                %Use previous weights as initial conditions to improve speed
                w0 = nanmean(weights([true oc],1:ii-1,kk),2);
                if ~any(isnan(w0))
                    w = ComputePeerPredictionRoot(q,[ones(size(pop,1),1) pop(:,oc)],dt,w0);
                else
                    w = ComputePeerPredictionRoot(q,[ones(size(pop,1),1) pop(:,oc)],dt);
                end
                
                temp = nan(nbC+1,1);
                
                temp([true oc]) = w;
                weights(:,ii,kk) = temp;
            end
            warning on
            
            
            %Now compute prediction
            ep  = crix(ii,:);
            
            bints = ts1(ep(1):ep(2));
            N = length(bints);
            popt = nan(N,nbC);
             %  qpopt = nan(N,nbC);
            for jj = 1:nbC
                 spk = cell2mat(Epoch2Mat(spikes.times{jj},ts1(ep)'));
              %  popt(:,jj) = nanconvn(histc(spk,bints),k);
                 popt(:,jj) = histc(spk,bints);
            end
            
            
            qt = popt(:,kk);
            r = sum(qt)/length(qt)/dt;
            LogFun = @(x) SpkTrainLogLikelihood(qt,  modifiedExp([ones(size(popt,1),1) popt(:,oc)]*x),dt) - 0.25*x'*x;
            Lf = LogFun(w);
            L0 = SpkTrainLogLikelihood(qt,repmat(r,[size(popt(:,oc),1) 1]),dt);
            
            %             disp(Lf-L0)
            L(ii,kk) = Lf-L0;
            
            
        end
    end
end

fprintf('done\n')
totlen = length(ts1)*dt;


L = nansum(L);
L = L/totlen;
L = L(:);

weights = squeeze(nanmean(weights,2));

end