function CM = graph_constructionKLS(A, method, parcellation)
% KL divergence method following Kong et al. 2014
%
% INPUT: nii images (already smoothed, it's another preprocessing step of
% course which will always be done first) 
%
% Steps: 
%   1. nk_WriteVol --> save image to disk (.nii); compute network; delete
%       this temp.nii; loop through images 
%   within Python script
%       1. compute probability density functions of the gray matter volumes 
%           for each region of the provided atlas 
%       2. compute symmetric KL divergence between each pair of probability
%           density functions
%       3. bring data in the right format (flat connectivity matrix, upper
%           triangle); dimensions: n people x (n/2-n) edges; either safe
%           this vector to csv-file or add to existing one 
% 
% OUTPUT: matrix format, dimensions: n people x (n/2-n) edges 
% 
% --> adapted Python script from my thesis
% potential issue: I also used R scripts for the computation ...





end