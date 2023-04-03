function [pY, mapping] = cv_PerfICA(Y, varargin)
% function to extract independent components from data 
% Input: Y - input data (N*D) 
% Output: pY - output data (N*IC) 
% in training mode, the ICs are computed; during testing the ICs are
% projected onto the unseen data (after unseen data is whitened)
global MODELDIR

% in training mode 
if nargin <= 1

    pY = pyrunfile('cv_py_PerfICA.py', 'ica_model', 'S', ...
        mode = 'train', ...
        data = Y, ...
        num_ics = 10, ...
        rootdir = MODELDIR); 
  

    pY = icasig'; % transpose so that nrows = nsamples, ncols = nICs
    mapping.ica_model = ica_model;
    
else % test mode
    mapping = varargin{1};
    pY = pyrunfile('cv_py_PerfICA.py', 'S' , ...
        mode = 'test', ...
        ica_model = mapping.ica_model, ...
        data = Y); 
 
end


