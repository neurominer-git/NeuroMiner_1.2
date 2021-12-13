% ==========================================================================
% FORMAT [param, model] = nk_GetParam_MEXELM(Y, label, SlackParam, ~, ...
%                                           ModelOnly)
% ==========================================================================
% Train LIBLINEAR models and evaluate their performance using Y & label, 
% SlackParam,
% if ModelOnly = 1, return only model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 11/2020

function [param, model] = nk_GetParam2_LIBLIN(Y, label, ModelOnly, cmd)
global SVM EVALFUNC CMDSTR MODEFL                          

param = [];

% Check if weighting is necessary
cmd = nk_SetWeightStr(SVM, MODEFL, CMDSTR, label, cmd);

if iscell(Y) 
    % MKL-based learning not implemented yet
   
else % Univariate case
    
    model = train_liblin22(label, sparse(Y), cmd);
	model = nk_BuildCalibrationModel(SVM, MODEFL, model, Y, label);
	
    %fprintf('\n%s',cmdstr)
    if ~ModelOnly
        [param.target, ...
            param.val, ...
            param.dec_values] = predict_liblin22(label, sparse(Y), model, SVM.LIBLIN.b);
        if SVM.RVMflag, param.dec_values = nk_CalibrateProbabilities(param.dec_values); end
        param.val = EVALFUNC(label, param.dec_values);
    end
end