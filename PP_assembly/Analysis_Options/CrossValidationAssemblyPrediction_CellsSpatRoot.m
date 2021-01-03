function [ses] = CrossValidationAssemblyPrediction_CellsSpatRoot(root,tempSmoothing,compPeers)

%
% [L,weights] = CrossValidationAssemblyPrediction_Cells(Starget,Qassembly,Qtarget,ep)
%
% computs the loglikelihood prediction of the spike train Starget from the
% binned spike train matrix of its peers (Qassembly). The prediction is built
% upon a Generalized Linear Model of the multivariate data given by
% Qassembly. The algorithm searches for a set of weights that maximize the
% prediction of the binned spike train Qtarget (the training set) and then
% computes the likelihood (during the test set). Everything is the same as
% in Harris et al., Nature, 2004.
%

if ~exist('compPeers','var') || isempty(compPeers)
  compPeers = 1;
end

% Parameters:
nbEp = length(root.ts); % number of epochs (i.e., cross valitation chunks); Epoch should be pre-set, so the number of cells in ts is equal to the number of cross-validations

% prep linear position
t = root.user_def.tracking;
root.b_myvar = mod(t.theta,2*pi);

% prep theta
if isempty(root.active_lfp)
  % theta phase will be determined by the LFP recorded on the tetrode with
  % the largest number of cells as the likely pyramidal layer theta
  lfpInd = mode(root.cel(:,1));
  root = root.LoadLFP(lfpInd,'downsample',250);
  root.active_lfp = lfpInd;
  if isempty(root.active_lfp), error('LFP failed to load'); end
  root.b_lfp(root.active_lfp) = root.b_lfp(root.active_lfp).AppendThetaPhase;
else
  error('This needs to be coded')
end

% bins = 0:1:130; % Sam's code. Spatial bins to split linearized position into
nBins = 345; % spatial bins
bins = linspace(0,2*pi,nBins); % ELN circle track is about 345 cm long
cel = root.cells;
l = 0.05; % smoothing width for smoothing ratmap & phase modulation estimates
dt = .0032;               % temporal sampling rate taken from Harris et al

fprintf('Launching Cross-validated Peer info\n')

% initializations
nbC = size(root.cells,1); % number of cells
Lf_trn = NaN(nbEp,nbC);  % Likelihood score, on training data, for [population + spatial]
La_trn = NaN(nbEp,nbC);  % Likelihood score, on training data, for [spatial + theta angle]
Ls_trn = NaN(nbEp,nbC);  % Likelihood score, on training data, for [spatial]
L0_trn = NaN(nbEp,nbC);  % Likelihood score, on training data, for [mean firing rate]
Lf_tst = NaN(nbEp,nbC);  % Likelihood score, on testing data, for [population + spatial]
La_tst = NaN(nbEp,nbC);  % Likelihood score, on testing data, for [spatial + theta angle]
Ls_tst = NaN(nbEp,nbC);  % Likelihood score, on testing data, for [spatial]
L0_tst = NaN(nbEp,nbC);  % Likelihood score, on testing data, for [mean firing rate]
L_fVs_tst = NaN(nbEp,nbC); % gain with position & peers over just spatial position on testing data
L_fVa_tst = NaN(nbEp,nbC); % gain with population & peers over spatial & theta angles [NOTE: not a reasonable comparison]
L_aVs_tst = NaN(nbEp,nbC); % gain with positing and theta angle over spatial position on testing data
L_sV0_tst = NaN(nbEp,nbC); % gain with spatial position over avg firing rate on testing data on testing data
L_aV0_tst = NaN(nbEp,nbC); % gain with spatial postion & theta angle over average firing rate on testing data
L_fV0_tst = NaN(nbEp,nbC); % gain with peers & spatial over agerage firing rate on testing data


% temporal smoothing terms
if ~exist('tempSmoothing','var') || isempty(tempSmoothing)
  tempSmoothing = 0.020; % 20ms
end
tSigma = tempSmoothing / dt;
tKernel = gaussian2Dfilter([ceil(8*tSigma)+1 1],[tSigma 1]);
tKernel = tKernel ./ max(tKernel); % set max to 1 to preserve spike counts

weights = nan(nbC,nbEp,nbC);
eps = root.epoch;
for ii=1:nbEp
  fprintf('.')
  
  tr = ~ismember(1:nbEp,ii);
  root.epoch = eps(tr,:);
  
  % define time bins at new temporal scale
  bints = cellfun(@(a) a(1):dt:a(end),root.ts,'uni',0);
  
  % initializations
  N = length(cell2mat(bints')); % number of time bins to predict
  pop = nan(N,nbC); % number of spikes in each of N time bins for eacn of the nbC cells
  f = zeros(N,nbC); % predicted spiking of each cell based on place for each time bin
  ft = zeros(N,nbC); % predicted spiking of each cell based on ratemap and theta modulation for each time bin
  
  % compute ratemaps for each cell
  [rm,vm_theta,vm_kappa] = deal(nan(size(root.cel,1),length(bins)));
  loc = cellfun(@(a,b,c) interp1(a,b,c),root.ts,root.myvar,bints,'uni',0); % location at new timescale

  for ind = 1:size(root.cel,1)
    spk_linpos = root.spk.myvar(:,ind); % eln code    
    [ratemap,spkcnt(ind,:),occ(ind,:)] = ratemap_Harris2003style(cat(2,loc{:})',cat(1,spk_linpos{:}),bins,l);
    rm(ind,:) = ratemap/dt;
    
    % compute phase locking prefs for each cell at each position w von Mises distro
    spkTh = root.spk.theta(:,ind);
    [vm_theta(ind,:),vm_kappa(ind,:)] = vonMises_Harris2003style(cat(1,spkTh{:}),cat(1,spk_linpos{:}),bins,l);
  end
  
  
  
  % position of animal in terms of bin number of discretized linearized position at new temporal scale
  pos = cellfun(@(a) cumsum([a(1) circDiff(a)]), root.myvar, 'uni',0); % position must be unwrapped to avoid interpolation errors around the start point of the circle track
  loc = cellfun(@(a,b,c) interp1(a,b,c),root.ts,pos,bints,'uni',0); % resample timebins of position
  [~,loc] = cellfun(@(a) histc(mod(a,2*pi),bins),loc,'uni',0); % match positions to bins of ratemap
  loc = cell2mat(loc');
  % set out of bounds positions to first bin
  badloc = loc==0;
  loc(badloc) = 1;
  
  
  % theta phase at new temporal scale
  thPhs = cellfun(@(a) cumsum([a(1) circDiff(a)]), root.lfp.theta_phase, 'uni',0); % theta phase must be unwrapped to avoid interpolation errors around 2pi
  thPhs = cellfun(@(a,b,c) interp1(a,b,c), root.lfp.ts, thPhs, bints, 'uni',0);    % match to new temporal scale
  thPhs = mod(cell2mat(thPhs'),2*pi); % concatonate and rewrap to circle
  
  
  %get poulation
  for jj = 1:nbC
    % get the spike counts per time bin for each cell
    spk = root.spk_ts(root.cells(jj,:));
    nspE = cellfun(@(a) numel(a),spk); spk(nspE==0) = {-1}; % eln code
    spk = cellfun(@(a,b) linearize(histc(a,b)),spk,bints,'uni',0); % this will fail with a size mismatch error if there are not spikes in every epoch bin
    spk = cellfun(@(a) nanconvn(a,tKernel,'1d'),spk,'uni',0); % smooth population activity
    pop(:,jj) = cell2mat(spk);
    
    % predict firing of each cell based on ratemap given the trajectory
    temp = rm(jj,loc);
    temp(badloc) = 0;
    f(:,jj) = temp;
    
    % predict firing of each cell based on theta phase & ratemap given the trajectory
    K = vm_kappa(jj,loc); 
    temp = exp(K.*cos(thPhs-vm_theta(jj,loc)))./besseli(0,K); %equation 8 from Harris supplementals
    temp(badloc) = 0;
    ft(:,jj) = f(:,jj).*temp';
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Now compute predictions %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%
  root.epoch = eps(ii,:);
  
  bints = root.ts(1):dt:root.ts(end);
  N = length(bints);
  f_t = zeros(N,nbC);
  ft_t = zeros(N,nbC);
  pop_t = nan(N,nbC);
  
  % compute position at new temporal scale
  pos_t = cumsum([root.myvar(1) circDiff(root.myvar)]); % position must be unwrapped to avoid interpolation errors around the start point of the circle track
  loc_t = interp1(root.ts,pos_t,bints); % resample timebins of position
  [~,loc_t] = histc(mod(loc_t,2*pi),bins); % rewrap to circle then match to bins of ratemap
  badloc_t = loc_t==0;
  loc_t(badloc_t) = 1;
  
  % compute theta phase at new temporal scale
  thPhs_t = root.lfp.theta_phase;
  thPhs_t = cumsum([thPhs_t(1) circDiff(thPhs_t)]); % theta phase must be unwrapped to avoid interpolation errors around 2pi
  thPhs_t = interp1(root.lfp.ts,thPhs_t,bints);    % match to new temporal scale
  thPhs_t = mod(thPhs_t,2*pi); % concatonate and rewrap to circle
  
  for jj = 1:nbC
    % number of actual spikes per time bin
    spk = root.spk_ts(root.cells(jj,:)); if isempty(spk), spk = -1; end
    pop_t(:,jj) = histc(spk,bints);
    
    % predicted firing rates based on position for each cell
    temp = rm(jj,loc_t);
    temp(badloc_t) = 0;
    f_t(:,jj) = temp;
    
    % predict firing of each cell based on theta phase & ratemap given the trajectory
    temp = circ_vmpdf(thPhs_t,vm_theta(jj,loc_t),vm_kappa(jj,loc_t));
    temp(badloc_t) = 0;
    ft_t(:,jj) = f_t(:,jj).*temp;
  end
  
  for kk = 1:nbC
    q = pop(~badloc,kk); % target cell spike counts from training set
    oc = ~ismember(1:nbC,kk); % other cell's indices
    if any(q>0)
      
      % perform fitting on training set
      warning off
      if ii==1
        temp = nan(nbC,1);
        w = ComputePeerPredictionSpatRoot(q,pop(~badloc,oc),f(~badloc,kk),dt);
        
        temp(oc) = w;
        weights(:,ii,kk) = temp;
      else
        %Use previous weights as initial conditions to improve speed
        w0 = nanmean(weights(oc,1:ii-1,kk),2);
        if ~any(isnan(w0))
          w = ComputePeerPredictionSpatRoot(q,pop(~badloc,oc),f(~badloc,kk),dt,w0);
        else
          w = ComputePeerPredictionSpatRoot(q,pop(~badloc,oc),f(~badloc,kk),dt);
        end
        
        temp = nan(nbC,1);
        
        temp(oc) = w;
        weights(:,ii,kk) = temp;
      end
      warning on
      
      % check performance on 'training set'
      % predict based on population and position
      Lf_trn(ii,kk) = SpkTrainLogLikelihood(q,dt*modifiedExp(pop(~badloc,oc)*w).*f(~badloc,kk))/(length(q)*dt);
      % predict based on position and theta phase
      La_trn(ii,kk) = SpkTrainLogLikelihood(q,dt*ft(~badloc,kk))/(length(q)*dt); % likelihood given angle of theta
      % predict based on position
      Ls_trn(ii,kk) = SpkTrainLogLikelihood(q,dt*f(~badloc,kk))/(length(q)*dt);
      % null model based on average firing rate alone
      L0_trn(ii,kk) = SpkTrainLogLikelihood(q,repmat(mean(q),size(q)))/(length(q)*dt);     
      
      % check performance on testing set
      qt = pop_t(~badloc_t,kk); % spiking of current cell
      % predict based on population and position
      Lf_tst(ii,kk) = SpkTrainLogLikelihood(qt,dt*modifiedExp(pop_t(~badloc_t,oc)*w).*f_t(~badloc_t,kk))/(length(qt)*dt);
      % predict based on position and theta phase
      La_tst(ii,kk) = SpkTrainLogLikelihood(qt,dt*ft_t(~badloc_t,kk))/(length(qt)*dt); % likelihood given angle of theta
      % predict based on position
      Ls_tst(ii,kk) = SpkTrainLogLikelihood(qt,dt*f_t(~badloc_t,kk))/(length(qt)*dt);
      % null model based on average firing rate alone
      L0_tst(ii,kk) = SpkTrainLogLikelihood(qt,repmat(mean(q),size(qt)))/(length(qt)*dt);
      
      % information rates per spike
      L_fVs_tst(ii,kk) = (Lf_tst(ii,kk)-Ls_tst(ii,kk)); % gain with position & peers over just position
      L_fVa_tst(ii,kk) = (Lf_tst(ii,kk)-La_tst(ii,kk)); % gain with position & peers over position & theta angle
      L_aVs_tst(ii,kk) = (La_tst(ii,kk)-Ls_tst(ii,kk)); % gain with theta phase & position over avg firing rate
      L_sV0_tst(ii,kk) = (Ls_tst(ii,kk)-L0_tst(ii,kk)); % gain with position over avg firing rate
      L_aV0_tst(ii,kk) = (La_tst(ii,kk)-L0_tst(ii,kk)); % gain with position & theta angle over avg firing rate
      L_fV0_tst(ii,kk) = (Lf_tst(ii,kk)-L0_tst(ii,kk)); % gain with position & population over avg firing rate
      
      nSpks_trn(ii,kk) = sum(q);
      nSpks_tst(ii,kk) = sum(qt);
    end
  end
end

fprintf('done\n')


% Colapse across n-1
Lf_trn = nanmean(Lf_trn);
La_trn = nanmean(La_trn);
Ls_trn = nanmean(Ls_trn);
L0_trn = nanmean(L0_trn);
Lf_tst = nanmean(Lf_tst);
La_tst = nanmean(La_tst);
Ls_tst = nanmean(Ls_tst);
L0_tst = nanmean(L0_tst);
L_fVs_tst = nanmean(L_fVs_tst);
L_fVa_tst = nanmean(L_fVa_tst);
L_aVs_tst = nanmean(L_aVs_tst);
L_sV0_tst = nanmean(L_sV0_tst);
L_aV0_tst = nanmean(L_aV0_tst);
L_fV0_tst = nanmean(L_fV0_tst);

weights = squeeze(nanmean(weights,2));
nSpks_trn = nanmean(q);
nSpks_tst = nanmean(qt);

% Package it up to easy shipping
ses.Lf_trn = Lf_trn(:);
ses.La_trn = La_trn(:);
ses.Ls_trn = Ls_trn(:);
ses.L0_trn = L0_trn(:);
ses.Lf_tst = Lf_tst(:);
ses.La_tst = La_tst(:);
ses.Ls_tst = Ls_tst(:);
ses.L0_tst = L0_tst(:);
ses.L_fVs_tst = L_fVs_tst(:);
ses.L_fVa_tst = L_fVa_tst(:);
ses.L_aVs_tst = L_aVs_tst(:);
ses.L_sV0_tst = L_sV0_tst(:);
ses.L_aV0_tst = L_aV0_tst(:);
ses.L_fV0_tst = L_fV0_tst(:);

ses.weights = weights;
ses.cel = cel;
ses.nSpks_trn = nSpks_trn;
ses.nSpks_tst = nSpks_tst;
ses.rm = rm;
ses.vm_theta = vm_theta;
ses.vm_kappa = vm_kappa;

