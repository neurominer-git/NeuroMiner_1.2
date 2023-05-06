function [pY, mapping] = cv_PerfICA(Y, opt, mode)
% function to extract independent components from data 
% Input: Y - input data (N*D) 
% Output: pY - output data (N*IC) 
% in training mode, the ICs are computed; during testing the ICs are
% projected onto the unseen data (after unseen data is whitened)
global MODELDIR

if strcmp(mode, 'train') % in training mode 

    if isfield(opt, 'dims')
        n_ics = opt.dims; % potentially add more options here, then opt should be struct with fields
    else 
        n_ics = 0;
    end


    pY = pyrunfile('cv_py_PerfICA.py', 'ica_model', 'S', ...
        mode = mode, ...
        data = Y, ...
        num_ics = n_ics, ...
        rootdir = MODELDIR); 
  

    pY = S'; % transpose so that nrows = nsamples, ncols = nICs
    mapping.ica_model = ica_model;
    
else % test mode
    mapping = varargin{1};
    pY = pyrunfile('cv_py_PerfICA.py', 'S' , ...
        mode = mode, ...
        ica_model = mapping.ica_model, ...
        data = Y); 
 
end


