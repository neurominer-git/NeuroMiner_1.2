% =========================================================================
% FORMAT param = nk_MLI_config(param)
% =========================================================================
% NeuroMiner menu configurator for the interpretation of model predictions 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 05 / 2022
function [param, act] = nk_MLI_config(param)
global NM

upper_thresh = 95;
lower_thresh = 5;
nperms  = 1000;
max_iter = 1000;
n_visited = 100;
frac    = .1;
method  = 'posneg';
if isfield(param,'MULTI') && param.MULTI.flag
    multiflag = true; groupmode = 2; 
else
    multiflag = false;
    groupmode = 1; 
end
trainwithCV2Ts = 1;

if ~isfield(NM.TrainParam,'MLI')
    NM.TrainParam.MLI.method = method;
    NM.TrainParam.MLI.nperms = nperms;
    NM.TrainParam.MLI.frac = frac;
    NM.TrainParam.MLI.upper_thresh = upper_thresh;
    NM.TrainParam.MLI.lower_thresh = lower_thresh;
    NM.TrainParam.MLI.max_iter = max_iter;
    NM.TrainParam.MLI.n_visited = n_visited;
end

switch method
    case 'posneg'
        MethodStr = sprintf('Upper/lower percentiles (%g/%g)', upper_thresh, lower_thresh);
        OcclusionUpperThreshStr = sprintf('Define upper percentile [%g]|', upper_thresh);
        OcclusionLowerThreshStr = sprintf('Define lower percentile [%g]|', lower_thresh);

    case 'median'
        MethodStr = 'Median';
        OcclusionUpperThreshStr = '';
        OcclusionLowerThreshStr = '';
end
OcclusionMethodStr = ['Define ncclusion method [ ' MethodStr ' ]|'];

mnustr = [OcclusionMethodStr ...
          OcclusionUpperThreshStr ...
          OcclusionLowerThreshStr ...
          ];

act = nk_input('MAIN INTERFACE >> DEFINE PARAMETERS >> Model interpreter settings',0,'mq', ...
                [ ...
                 '# of iterations [ ' o.npermsstr ' ]|' ...
                 'Model retraining mode [ ' o.trainwithCV2Tsstr ' ]'],1:3);

nk_PrintLogo


switch act

    case 1
        meanflag = nk_input('Ensemble construction mode', 0, 'm', ...
                                ['Aggregate ALL base learners into ONE ensemble|' ...
                                 'Compute mean decisions of CV1 partion ensembles before aggregating'], ...
                                 [1,2], meanflag);
    case 2
        groupmode = nk_input('Define group processing mode',0, 'm', ...
                                ['OOCV prediction at binary predictors'' optimum parameters|' ...
                                 'OOCV prediction at multi-group predictors'' optimum parameters|' ...
                                 'OOCV prediction at binary & multi-group predictors'' optimum parameters'], ...
                                 1:3, groupmode);
            
    case 3
        trainwithCV2Ts = nk_input('Model retraining mode', 0, 'm', ...
                                 ['As defined for CV training phase|' ...
                                  'Use all available CV data (CV1+CV2)'], ...
                                  [1,2], trainwithCV2Ts); 
end

param.OOCV.meanflag         = meanflag;
param.OOCV.groupmode        = groupmode;
param.OOCV.trainwithCV2Ts   = trainwithCV2Ts;
param.OOCV.savemodels       = savemodels;
param.OOCV.saveoocvdata     = saveoocvdata;

end