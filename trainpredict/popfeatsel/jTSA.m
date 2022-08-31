%[2015]-"TSA: Tree-seed algorithm for continuous optimization"

% (9/12/2020)

function [Sf, curve] = jTSA(rx, Ps, N, max_Iter, lb, ub, thres, ST)

% Objective function
fun = @jFitness; 
funoper = str2func(rx.evaldir2); 
funopti = str2func(rx.minmaxstr);

% Number of dimensions
dim = size(rx.Y,2); 

% Initial (5)
X   = zeros(N,dim); 
for i = 1:N
  for d = 1:dim
    X(i,d) = lb + (ub - lb) * rand();
  end 
end
% Fitness
fit = zeros(1,N); 
for i = 1:N 
  fit(i) = fun(rx, X(i,:) > thres , Ps);
end
% Best solution (6)
[fitG, idx] = funopti(fit);
Xgb         = X(idx,:);

% Maximum & minimum number of seed
Smax = round(0.25 * N); 
Smin = round(0.1 * N);

% Pre
curve = zeros(1,max_Iter); curve(1) = fitG;
t = 2;
% Iteration
while t <= max_Iter
  for i = 1:N
    % Random number of seed
    num_seed = round(Smin + rand()* (Smax - Smin)); 
    Xnew     = zeros(num_seed, dim);
    for j = 1:num_seed
      % Random select a tree, but not i
      RN = randperm(N); 
      RN(RN == i) = []; 
      r  = RN(1);
      for d = 1:dim
        % Alpha in [-1,1]
        alpha = -1 + 2 * rand();
        if rand() < ST  
          % Generate seed (3)
          Xnew(j,d) = X(i,d) + alpha * (Xgb(d) - X(r,d));
        else
          % Generate seed (4)
          Xnew(j,d) = X(i,d) + alpha * (X(i,d) - X(r,d));
        end
      end
      % Boundary
      XB = Xnew(j,:); XB(XB > ub) = ub; XB(XB < lb) = lb;
      Xnew(j,:) = XB;
    end
    % Fitness
    for j = 1:num_seed
      % Fitness
      Fnew = fun( rx , Xnew(j,:) > thres , Ps);
      % Greedy selection
      if funoper(Fnew, fit(i))
        fit(i) = Fnew;
        X(i,:) = Xnew(j,:);
      end
    end
  end
  % Best solution (6)
  [fitG_new, idx] = funopti(fit);
  Xgb_new         = X(idx,:);
  % Best update
  if funoper(fitG_new, fitG)
    fitG = fitG_new;
    Xgb  = Xgb_new;
  end
  % Store 
  curve(t) = fitG;
  %fprintf('\nIteration %d Best (TSA)= %f',t,curve(t))
  t = t + 1;
end
% Select features
Pos   = 1:dim; 
Sf    = Pos((Xgb > thres) == 1); 