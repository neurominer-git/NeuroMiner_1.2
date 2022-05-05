function [param, model] = nk_GetParam2_RNDFOR(Y, label, ModelOnly, Param)
global EVALFUNC 

param = [];
       
if iscell(Y) 
    % MKL-based learning not implemented yet
   
else % Univariate case
    
    model = pyrunfile('nm_RF_Classifier_train.py', 'model_file', ...
        Y, label, int64(Param(1)), int64(Param(2)));

    %fprintf('\n%s',cmdstr)
%     if ~ModelOnly
%         [param.target] = predict_liblin(label, Y, model);
%         param.dec_values = param.target; 
%         param.val = EVALFUNC(label, param.dec_values);
%     end
end
