% =========================================================================
% FORMAT param = AUC(expected, predicted)
% =========================================================================
% Compute Area-Under-the-Curve for binary classification problems
% Based on the code of Chih-Jen Lin
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 08/2015

function param = AUC(expected, predicted)
param = []; 
if isempty(expected), return; end
ind0 = expected ~=0;
expected = expected(ind0); 
predicted = predicted(ind0);
uE = unique(expected);
if numel(uE)==2
      param = fastAUC(expected, predicted,1);
else
    error('AUC works only for binary problems'); 
end