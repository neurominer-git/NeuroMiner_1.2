function [Sf,curve] = jPSO(r, Ps , opts)
% -------------------------------------------------------------------------
% FORMAT function [Sf,curve] = jPSO(r, Ps , opts)
% -------------------------------------------------------------------------
% Code of the particle swarm optimization algorithm downloaded from:
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
% Refs
% [1995]-"Particle Swarm Optimization" 
% [1998]-"A modified particle swarm optimizer"
% =========================================================================
% (c) Nikolaos Koutsouleris, 03/2022

% prepare
if ~exist('opts','var') || isempty(opts)
    opts = nk_PSO_config([],1);
end
N = opts.N;
max_Iter = opts.max_Iter;
c1 = opts.c1;
c2 = opts.c2;
w = opts.w;
lb = opts.lb;
ub = opts.ub;
thres = opts.thres;
Vmax = (ub-lb)/2; % Maximum velocity

% Objective function
fun = @jFitness; 
funop = str2func(r.evaldir2); 
% Number of dimensions
dim = size(r.Y,2); 
% Initial 
X   = zeros(N,dim);
V   = zeros(N,dim); 
for i = 1:N
  for d = 1:dim
    X(i,d) = lb + (ub - lb) * rand();
  end
end  
% Fitness
fit  = zeros(1,N); 
fitG = r.optparam;
for i = 1:N 
  fit(i) = fun(r, X(i,:) > thres, Ps); 
  % Gbest update
  if funop(fit(i) , fitG)
    Xgb  = X(i,:); 
    fitG = fit(i);
  end
end
% PBest
Xpb  = X; 
fitP = fit;
% Pre
curve = inf;
t = 1;  
% Iterations
while t <= max_Iter
  for i = 1:N
    for d = 1:dim
      r1 = rand();
      r2 = rand();
      % Velocity update 
      VB = w * V(i,d) + c1 * r1 * (Xpb(i,d) - X(i,d)) + ...
        c2 * r2 * (Xgb(d) - X(i,d)); 
      % Velocity limit
      VB(VB > Vmax) = Vmax;  VB(VB < -Vmax) = -Vmax;
      V(i,d) = VB;
      % Position update 
      X(i,d) = X(i,d) + V(i,d);
    end
    % Boundary
    XB = X(i,:); XB(XB > ub) = ub; XB(XB < lb) = lb; 
    X(i,:) = XB;
    % Fitness
    fit(i) = fun( r , X(i,:) > thres, Ps );
    % Pbest update
    if funop(fit(i), fitP(i)) 
      Xpb(i,:) = X(i,:); 
      fitP(i)  = fit(i);
    end
    % Gbest update
    if funop(fitP(i), fitG)  
      Xgb  = Xpb(i,:);
      fitG = fitP(i);
    end
  end
  curve(t) = fitG;
  t = t + 1;
end

% Select features based on selected index
Pos   = 1:dim;
Sf    = Pos((Xgb > thres) == 1); 

