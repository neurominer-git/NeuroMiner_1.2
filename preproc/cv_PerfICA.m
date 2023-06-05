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

    [model_file, S, ICs] = pyrunfile('cv_py_PerfICA.py', ["model_file" "S", "ICs"], ...
        mode = 'train', ...
        data = Y, ...
        num_ics = int64(n_ics), ...
        rootdir = MODELDIR); 
  
    pY = py2mat(S); % transpose so that nrows = nsamples, ncols = nICs
    mapping.ica_py_model_file = model_file;
    mapping.vec = py2mat(ICs)';
    mapping.sampleMean = mean(data,1);
elseif strcmp(mode, 'test') % test mode
    mapping = opt;
    S = pyrunfile('cv_py_PerfICA.py', 'S' , ...
        mode = 'test', ...
        ica_model = mapping.ica_py_model_file, ...
        data = Y); 
    pY = py2mat(S);
elseif strcmp(mode, 'inverse_transform')
    mapping = opt; 
    S = pyrunfile('cv_py_PerfICA.py', 'S' , ...
        mode = 'test', ...
        ica_model = mapping.ica_py_model_file, ...
        data = Y); 
    pY = py2mat(S);
end


