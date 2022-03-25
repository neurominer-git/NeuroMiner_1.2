function [ Sf, curve ] = jPFA(rx, Ps, opts)
% -------------------------------------------------------------------------
% FORMAT function [Sf,curve] = jPFA(r, Ps , opts)
% -------------------------------------------------------------------------
% Code of the pathfinder algorithm downloaded from:
% https://github.com/JingweiToo/Wrapper-Feature-Selection-Toolbox/
% ... and modified for use in NeuroMiner 
%
% Inputs:
% r         : Data input structure as provided by rfe_algo_settings
% Ps        : wrapped (learning) algorithm's parameters
% opts      : wrapper optimization parameters
%
% Outputs:
% Sf        : Selected features (vector with featur indices)
% curve     : optimization history
%
% Refs:
% [2019]-"A new meta-heuristic optimizer: Pathfinder algorithm"
% https://www.sciencedirect.com/science/article/pii/S1568494619301309
% =========================================================================
% (c) Nikolaos Koutsouleris, 03/2022

% prepare
if ~exist('opts','var') || isempty(opts)
    opts = nk_PFA_config([],1);
end
N = opts.N;
max_Iter = opts.max_Iter;
lb = opts.lb;
ub = opts.ub;
thres = opts.thres;

% Objective function
fun = @jFitness; 
funop = str2func(rx.evaldir2); 

% Number of dimensions
dim = size(rx.Y,2); 

% Initial 
X   = zeros(N,dim); 
for i = 1:N
	for d = 1:dim
    X(i,d) = lb + (ub - lb) * rand();
	end
end
% Fitness
fit  = zeros(1,N); 
fitP = inf;
for i = 1:N
  fit(i) = fun(rx , X(i,:) > thres, Ps);
  % Pathfinder update
  if funop(fit(i), fitP)
    fitP = fit(i);
    Xpf  = X(i,:);
  end
end

% Set previous pathfiner
Xpf_old = Xpf; 

% Pre
Xpf_new = zeros(1,dim);
Xnew    = zeros(N,dim);

curve = zeros(1,max_Iter);
curve(1) = fitP;
t = 2;

% Iterations
while t <= max_Iter
  % Alpha & beta in [1,2]
  alpha = 1 + rand();
  beta  = 1 + rand();
  for d = 1:dim  
    % Define u2 in [-1,1]
    u2 = -1 + 2 * rand();
    % Compute A (2.6)
    A  = u2 * exp(-(2 * t) / max_Iter);
    % Update pathfinder (2.4) 
    r3 = rand();
    Xpf_new(d) = Xpf(d) + 2 * r3 * (Xpf(d) - Xpf_old(d)) + A;
  end
  % Boundary
  Xpf_new(Xpf_new > ub) = ub; Xpf_new(Xpf_new <lb) = lb;
  % Update previous path
  Xpf_old = Xpf;
  % Fitness
  Fnew = fun(rx, Xpf_new > thres, Ps);
  % Greedy selection
  if funop(Fnew, fitP)
    fitP = Fnew; 
    Xpf  = Xpf_new;
  end

  % Sort member
  [fit, idx] = sort(fit,rx.optfunc);
  X          = X(idx,:); 
  % Update first solution 
  if funop(Fnew, fit(1))
    fit(1) = Fnew; 
    X(1,:) = Xpf_new;
  end

  % Update 
  for i = 2:N
    % Distance (2.5)
    Dij = norm(X(i,:) - X(i-1,:));
    for d = 1:dim
      % Define u1 in [-1,1]
      u1  = -1 + 2 * rand();  
      % Compute epsilon (2.5)
      eps = (1 - (t / max_Iter)) * u1 * Dij;
      % Define R1, R2
      r1  = rand(); 
      r2  = rand(); 
      R1  = alpha * r1; 
      R2  = beta * r2;
      % Update member (2.3)
      Xnew(i,d) = X(i,d) + R1 * (X(i-1, d) - X(i,d)) + ...
        R2 * (Xpf(d) - X(i,d)) + eps;
    end
    % Boundary
    XB = Xnew(i,:); XB(XB > ub) = ub; XB(XB < lb) = lb;
    Xnew(i,:) = XB;
  end
  % Fitness
  for i = 2:N
    % Fitness
    Fnew = fun(rx, Xnew(i,:) > thres, Ps);
    % Selection
    if funop(Fnew, fit(i))
      fit(i) = Fnew; 
      X(i,:) = Xnew(i,:);
    end
    % Pathfinder update
    if funop(fit(i), fitP)
      fitP = fit(i);
      Xpf  = X(i,:);
    end
  end
  curve(t) = fitP;
  %fprintf('\nIteration %d Best (PFA)= %f',t,curve(t))
  t = t + 1;
end
% Select features
Pos   = 1:dim;
Sf    = Pos((Xpf > thres) == 1); 
