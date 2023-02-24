function [pY, IN.mpp] = cv_PerfICA(Y, varargin)
% function to extract independent components from data 
% Input: Y - input data (N*D) 
% Output: pY - output data (N*IC) 
% in training mode, the ICs are computed; during testing the ICs are
% projected onto the unseen data

% in training mode 
if nargin > 1
    [icasig, A, W] = fastica(Y');
    % Y = projection of original data on the ICs
    % A = ICs are represented in the columns 
    % W = separating matrix W can be used to project a new data on ICs (Yn =
    % W*Xn) 

    pY = icasig';
    IN.mpp.A = A; 
    IN.mpp.W = W;
else % test mode
    mapping = varargin(1);
    pY = mapping.W * Y';
end


