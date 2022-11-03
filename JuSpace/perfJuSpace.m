function [sY, IN] = perfJuSpace(Y, IN)
% =========================================================================

% =========================== WRAPPER FUNCTION ============================ 
    if ~exist('IN','var'), IN = []; end
    if iscell(Y) 
        sY = cell(1,numel(Y)); 
        for i=1:numel(Y), [sY{i}, IN] =  PerfJuSpace(Y{i}, IN); end
    else
        [ sY, IN ] = PerfJuSpace(Y, IN);
    end
end

% =========================================================================
function [Y, IN] = PerfJuSpace(Y, IN)
    
    inputStr = '';
    
    if ~isempty(IN.brainmask) && ~isempty(IN.atlas) && ~isempty(IN.cortype) && ~isempty(IN.autocorcorrect) && ~isempty(IN.petList)
        Y = apply_JuSpace(Y, IN.brainmask, IN.atlas, IN.cortype, IN.autocorcorrect, IN.petList);
    end
        
    
end
