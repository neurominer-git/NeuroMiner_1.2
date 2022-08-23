function [sY, IN] = graph_PerfGraphConstruction(Y, IN)
% =========================================================================

% =========================== WRAPPER FUNCTION ============================ 
    if ~exist('IN','var'), IN = []; end
    if iscell(Y) 
        sY = cell(1,numel(Y)); 
        for i=1:numel(Y), [sY{i}, IN] =  PerfGraphConstruction(Y{i}, IN); end
    else
        [ sY, IN ] = PerfGraphConstruction(Y, IN);
    end
end

% =========================================================================
function [Y, IN] = PerfGraphConstruction(Y, IN)

    %global VERBOSE SVM
% Defaults
%if isempty(IN),eIN=true; else eIN=false; end

%if eIN  

    if ~isempty(IN.method) 
        switch IN.method 
            case "KL divergence" 
                 if ~isempty(IN.parcellation) %&& (IN.parcellation == "Hammers.nii") 
                    R = graph_constructionKLS(Y, IN.method, IN.parcellation); % KLS = symmetric KL divergence method
                 end
            case "Normative network + 1"
                if ~isempty(IN.refGroup)
                    R = graph_constructionREFPLUSP(Y, IN.method, IN.refGroup, IN.simMeasure) %JBE = jackknife bias estimation method
                end
        end
    else
        R = Y;
    end
    
Y = R;
end