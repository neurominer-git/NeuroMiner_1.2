function [sY, IN] = perfROImeans(Y, IN)%, PREPROC, prevP)
% =========================================================================

% =========================== WRAPPER FUNCTION ============================ 
    if ~exist('IN','var'), IN = []; end
    if iscell(Y) 
        sY = cell(1,numel(Y)); 
        for i=1:numel(Y), [sY{i}, IN] =  PerfROImeans(Y{i}, IN)%, PREPROC, prevP); end
    else
        [ sY, IN ] = PerfROImeans(Y, IN)%, PREPROC, prevP);
    end
end

% =========================================================================
function [Y, IN] = PerfROImeans(Y, IN)%, PREPROC, prevP)
    
    inputStr = '';
    
    if isfield(IN, 'brainmask') && isfield(IN, 'atlas') 
        Y = compute_ROImeans(Y, IN.brainmask, IN.atlas)%, PREPROC, prevP);
    end
        
    
end