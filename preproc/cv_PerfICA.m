function [pY, mapping] = cv_PerfICA(Y, varargin)
% function to extract independent components from data 
% Input: Y - input data (N*D) 
% Output: pY - output data (N*IC) 
% in training mode, the ICs are computed; during testing the ICs are
% projected onto the unseen data (after unseen data is whitened)
global MODELDIR

% in training mode 
if nargin <= 1

    [model_file, S] = pyrunfile('cv_py_PerfICA.py', ["model_file" "S"], ...
        mode = 'train', ...
        data = Y, ...
        num_ics = 10, ...
        rootdir = MODELDIR); 
  

    pY = py2mat(S); % transpose so that nrows = nsamples, ncols = nICs
    mapping.ica_py_model_file = model_file;
    
else % test mode
    mapping = varargin{1};
    S = pyrunfile('cv_py_PerfICA.py', 'S' , ...
        mode = 'test', ...
        ica_model = mapping.ica_py_model_file, ...
        data = Y); 
    pY = py2mat(S);
 
end


