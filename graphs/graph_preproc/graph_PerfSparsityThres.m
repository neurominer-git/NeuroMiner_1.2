function [sY, IN] = graph_PerfSparsityThres(Y, IN)
% =========================================================================
% FORMAT function [sY, IN] = nk_PerfElimZeroObj(Y, IN)
% =========================================================================
% Remove features with zero-variance, and ANY Infs and NaNs.
% Furthermore, remove features with highly skewed distributions
% I/O arguments:
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 05/2018

% =========================== WRAPPER FUNCTION ============================ 
    if ~exist('IN','var'), IN = []; end
    if iscell(Y) 
        sY = cell(1,numel(Y)); 
        for i=1:numel(Y), [sY{i}, IN] =  PerfSparsityThres(Y{i}, IN); end
    else
        [ sY, IN ] = PerfSparsityThres(Y, IN );
    end
end

% =========================================================================
function [Y, IN] = PerfSparsityThres(Y, IN)

    %global VERBOSE SVM
% Defaults
%if isempty(IN),eIN=true; else eIN=false; end

%if eIN  

    if ~isempty(IN.perc)
        R = apply_sparsity_thres(Y, IN.p); 
    else
        R = Y;
    end
    
Y = R;
end