function P = nk_CalibrateProbabilities(P, cutoff)

if ~exist('cutoff','var') || isempty(cutoff)
    P = P - 0.5; 
else
    P = P - (cutoff+realmin);
end
