% ==========================================================================
% FORMAT [param, model] = nk_GetParam2_ELASVM(Y, Y, ModelOnly, Params)
% ==========================================================================
% Train Support vector elastic net model and evaluate their performance 
% if ModelOnly = 1, return only model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 04/2022
function [param, m] = nk_GetParam2_ELASVM( X, Y, ModelOnly, Params)
global EVALFUNC
% Parameters
param = [];
[~,p] = size(X); 
m = struct('t',Params(1),'lambda',Params(2),'classifier',Params(3),'e',Params(4));
m.C = 1/(2*m.lambda); 

%% Train LIBLINEAR model
% in high p low n setting use the primal formulation of LIBLINEAR
%if 2*p>n, m.classifier = 2; end
cmdstr = sprintf(' -c %5.10f -s %g -e %1.4f -q 1', m.C, m.classifier, m.e);

% Prepare new X and y
Xnew = [bsxfun(@minus,X,Y./m.t) bsxfun(@plus,X,Y./m.t)]';
Ynew = [ones(p,1);-ones(p,1)];

model = train_liblin244(Ynew, sparse(Xnew), cmdstr);
m.w = model.w; m.w = m.w';
m.alpha = m.C * max(1- Ynew.*(Xnew * m.w),0);
m.beta = m.t *(m.alpha(1:p) - m.alpha(p+1:2*p)) /sum(m.alpha);

if ~ModelOnly
    [param.target, param.dec_values] = nk_GetTestPerf_ELASVM([],X,[],m);
    param.val = EVALFUNC(label, param.target);
else
    param = [];
end
