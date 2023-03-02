function [pY, mapping] = cv_PerfICA(Y, varargin)
% function to extract independent components from data 
% Input: Y - input data (N*D) 
% Output: pY - output data (N*IC) 
% in training mode, the ICs are computed; during testing the ICs are
% projected onto the unseen data (after unseen data is whitened)

% in training mode 
if nargin <= 1

    % in the Matlab implementation of fastica(), whitening of the data is
    % internally performed, so we do not need to do it before hand 

    [S, A, W] = fastica(Y'); % default numOfIC = min(size(X)) TO DO: add numOfIC as hyperparameter to be specified in NM

    % S: This is a matrix of size ncomp x nsamples, 
    %   where ncomp is the number of independent components extracted and 
    %   nsamples is the number of samples in the input data. Each row of S 
    %   contains the time series of one independent component, which 
    %   represents a linear combination of the mixed signals. The 
    %   independent components are ordered in decreasing order of kurtosis.
    %
    % A: This is a matrix of size ncomp x nobs, where nobs is the number of 
    %   observed signals in the input data. Each row of A represents 
    %   the mixing matrix, which relates the observed signals to the 
    %   independent components. The rows of A are ordered in the same 
    %   way as the independent components in S.
    %
    % W: This is a matrix of size nobs x ncomp, which represents the 
    %   inverse of the mixing matrix. Each row of W corresponds to one 
    %   observed signal, and each column corresponds to one independent 
    %   component. The rows of W are ordered in the same way as the 
    %   rows of ica.A

    pY = S'; % transpose so that nrows = nsamples, ncols = nICs
    mapping.S = S;
    mapping.A = A; 
    mapping.W = W; % for testing context
else % test mode
    mapping = varargin{1};
    pY = mapping.W * Y'; % whiten test data
    pY = pY'; % transpose so that nrows = nsamples, ncols = nICs
end


