function [sY, IN] = perfCustomPreproc(Y, IN)
% =========================================================================

% =========================== WRAPPER FUNCTION ============================ 
    if ~exist('IN','var'), IN = []; end
    if iscell(Y) 
        sY = cell(1,numel(Y)); 
        for i=1:numel(Y), [sY{i}, IN] =  PerfCustomPreprocFunc(Y{i}, IN); end
    else
        [ sY, IN ] = PerfCustomPreprocFunc(Y, IN);
    end
end

% =========================================================================
function [Y, IN] = PerfCustomPreprocFunc(Y, IN)
    
    inputStr = '';
    
    if isfield(IN, 'ParamCombos')
        curParams = IN.ParamCombos(IN.p,:);
        Y = eval(sprintf("%s(Y, curParams);", IN.filename));
    else
        Y = eval(sprintf("%s(Y);",IN.filename));
    end
        
    
end