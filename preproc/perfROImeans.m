function [sY, IN] = perfROImeans(Y, IN)
% =========================================================================

% =========================== WRAPPER FUNCTION ============================ 
    if ~exist('IN','var'), IN = []; end
    if iscell(Y) 
        sY = cell(1,numel(Y)); 
        for i=1:numel(Y), [sY{i}, IN] =  PerfROImeans(Y{i}, IN); end
    else
        [ sY, IN ] = PerfROImeans(Y, IN);
    end
end

% =========================================================================
function [Y, IN] = PerfROImeans(Y, IN)
    
    inputStr = '';
    
    if isfield(IN, 'brainmask') && isfield(IN, 'atlas') 
        Y = compute_ROImeans(Y, IN.brainmask, IN.atlas);
    end
        
    
end