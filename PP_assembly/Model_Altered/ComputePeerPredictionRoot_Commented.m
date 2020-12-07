function [weights] = ComputePeerPredictionRoot_Commented(q,pop,dt,varargin)

% weights = ComputePeerPrediction(S,Q,ep)
% S: ts of the cell
% Q: binned spike train matrix (tsd object) of all other cells

% control weights, randomized
if ~isempty(varargin)
    weights0 = varargin{1};
else
    weights0 = randn(size(pop,2),1);
end

%fits population of binned spikes using a modified exponential function,
% and then compares to spike timestamps of a given cell?
L = @(x)-SpkTrainLogLikelihood(q,modifiedExp(pop*x),dt) + 0.25*x'*x; 
options = optimoptions('fminunc','display','off');

problem.objective = L;
problem.options = options;
problem.solver = 'fmincon';
problem.x0 = weights0;
problem.ObjectiveLimit = 1e-5;
problem.lb = -100*ones(1,size(pop,2));
problem.ub = 100*ones(1,size(pop,2));

weights = fmincon(problem); %fmincon = nonlinear program solver
%                          --so weights determined by log likelihood? 
% weights = fminunc(L,weights0)





