function [ optparam, optind , optfound, optmodel] = nk_WrapperFeatSelStrat(method, Y, label, Ynew, labelnew, Ps, FullFeat, FullParam, ActStr)
global TRAINFUNC VERBOSE RFE

optfound = false;
VERBOSE=true;
r = rfe_algo_settings(Y, label, Ynew, labelnew, Ps, FullFeat, FullParam, ActStr);

% Parameter setting
switch method
    case {'GA', 'PSO', 'PFA'}
        % PSO Algorithm
        [optind,curve] = feval(['j' method], r, Ps, RFE.Wrapper.(method)); 
    case 'TSA'
        % Parameters
        N        = 10; % number of seeds
        max_Iter = 100;
        lb       = 0;
        ub       = 1; 
        thres    = 0.5; 
        ST       = 0.1; 
        [optind,curve] = jTreeSeedAlgorithm(r, Ps, N, max_Iter, lb, ub, thres, ST); 

end
[~,optmodel] = TRAINFUNC(r.Y(:,optind), r.YL, 1, Ps); 
optparam = curve(end);

% Plot convergence curve
if VERBOSE
    fprintf(' %s: %g feats selected =>', method, numel(optind))
%     figure;
%     plot(1:max_Iter,curve);
%     xlabel('Number of generations');
%     ylabel('Fitness Value'); 
%     title(lb); grid on;
end
if feval(r.evaldir, optparam, r.optparam) 
     optfound = true; 
end