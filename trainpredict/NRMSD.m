% =========================================================================
% FORMAT param = NRMSD(expected, predicted)
% =========================================================================
% Compute Normalized Root of Mean Squared Deviation (NRMSD) of regression
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% nrmsd = 100 * [ rmse(sim, obs) / ( max(obs) - min(obs) ) ] 
% (c) Nikolaos Koutsouleris, 03/2022
function param = NRMSD(expected, predicted)
if isempty(expected), param = []; return; end
rmse = sqrt(MSE(expected,predicted));
param = rmse / range(expected) * 100;

end