function [ optparam, optind , optfound, optmodel] = GA(Y, label, Ynew, labelnew, Ps, FullFeat, FullParam, ActStr, N, max_Iter, CR, MR)
global TRAINFUNC VERBOSE
optfound = false;

r = rfe_algo_settings(Y, label, Ynew, labelnew, Ps, FullFeat, FullParam, ActStr);

% Parameter setting
if nargin<3
    N        = 10;
    max_Iter = 100;
    CR       = 0.8;
    MR       = 0.5;
end

% Genetic Algorithm
[~,optind,curve] = jGA1(r, Ps, N, max_Iter, CR, MR); 
[~,optmodel] = TRAINFUNC(r.Y(:,optind), r.YL, 1, Ps); 
optparam = curve(end);

% Plot convergence curve
if VERBOSE
figure;
    plot(1:max_Iter,curve);
    xlabel('Number of generations');
    ylabel('Fitness Value'); 
    title('GA'); grid on;
end
 if feval(r.evaldir, optparam, r.optparam) 
     optfound = true; 
 end