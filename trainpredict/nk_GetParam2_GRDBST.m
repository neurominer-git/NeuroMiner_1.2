% ==========================================================================
% FORMAT [param, model] = nk_GetParam_GRDBST(Y, label, ModelOnly, Params)
% ==========================================================================
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 03/2017

function [param, model] = nk_GetParam2_GRDBST(Y, label, ModelOnly, Params)
global EVALFUNC MODELDIR MODEFL

param = [];
options = nk_GenMatLearnOptions(Params);
%maxIters = uint32(options.maxIters);
%options.maxTreeDepth = uint32(options.maxTreeDepth);
%options.loss = 'squaredloss';
%model = SQBMatrixTrain( single(Y), label, maxIters, options );
switch MODEFL
    case 'classification'
        model = pyrunfile('py_classGRDBST_train.py', 'model_file', ...
            feat = Y, lab = label, ...
            n_est = int64(options.maxIters), ...
            l = string(options.loss), ...
            lr = double(options.learningRate), ...
            subsamp = double(options.subsamplingFactor),...
            n_maxdepth = int64(options.maxTreeDepth), ...
            rootdir = MODELDIR);
    case 'prediction'
        model = pyrunfile('py_regGRDBST_train.py', 'model_file', ...
            feat = Y, lab = label, ...
            n_est = int64(options.maxIters), ...
            l = string(options.loss), ...
            lr = double(options.learningRate), ...
            subsamp = double(options.subsamplingFactor),...
            n_maxfeat = int64(options.maxTreeDepth), ...
            rootdir = MODELDIR);
end

if ~ModelOnly
    [param.target, param.dec_values] = nk_GetTestPerf_GRDBST([], Y, label, model) ;
    param.val = EVALFUNC(label, param.target);
end
end