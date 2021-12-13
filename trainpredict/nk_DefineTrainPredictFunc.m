function [trainfunc, predictfunc] = nk_DefineTrainPredictFunc(experimental)

global SVM
% Create function strings
if exist('experimental','var') && ~isempty(experimental) && experimental
    TrainBase = 'nk_GetParam2_';
else
    TrainBase = 'nk_GetParam_';
end
PredictBase = 'nk_GetTestPerf_';
% Turn strings to function handles
trainfunc = str2func([ TrainBase SVM.prog ]);
predictfunc = str2func([ PredictBase SVM.prog ]) ;

end