% ==========================================================================
% FORMAT [param, model] = nk_GetParam2_LIBLIN(Y, label, ModelOnly, cmd)
% ==========================================================================
% Train LIBLINEAR models and evaluate their performance using Y & label, 
% SlackParam,
% if ModelOnly = 1, return only model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 04/2022

function [param, model] = nk_GetParam2_LIBLIN(Y, label, ModelOnly, cmd)
global SVM EVALFUNC CMDSTR MODEFL                          

param = [];

% Check if weighting is necessary
cmd = nk_SetWeightStr(SVM, MODEFL, CMDSTR, label, cmd);

if iscell(Y) 
    % MKL-based learning not implemented yet
   
else % Univariate case
    
    sY = sparse(Y);
    if isfield(SVM.LIBLIN,'AutoC') && SVM.LIBLIN.AutoC
        cmdopt = [' -C ' cmd ];
        OptC = train_liblin244(label, sY, cmdopt);
        c_old = extractBetween(cmd,"-c "," -");
        cmd = regexprep(cmd,c_old,sprintf("%g", OptC(1)));
    end
    model = train_liblin244(label, sY, cmd);
	model = nk_BuildCalibrationModel(SVM, MODEFL, model, Y, label);
	
    %fprintf('\n%s',cmdstr)
    if ~ModelOnly
        [param.target, ...
            param.val, ...
            param.dec_values] = predict_liblin244(label, sparse(Y), model, SVM.LIBLIN.b);
        if SVM.RVMflag, param.dec_values = nk_CalibrateProbabilities(param.dec_values); end
        param.val = EVALFUNC(label, param.dec_values);
    end
end