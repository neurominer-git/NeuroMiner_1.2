function CM = graph_construction(A, method, parcellation)
% KL divergence method following Kong et al. 2014
%
% INPUT: nii images (already smoothed, it's another preprocessing step of
% course which will always be done first) 
%
% Steps: 
%   1. compute probability density functions of the gray matter volumes 
%       for each region of the provided atlas 
%   2. compute symmetric KL divergence between each pair of probability
%       density functions
%   3. bring data in the right format (flat connectivity matrix, upper
%       triangle); dimensions: n people x (n/2-n) edges 
% 
% OUTPUT: matrix format, dimensions: n people x (n/2-n) edges 



end