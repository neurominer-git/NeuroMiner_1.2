function D = nk_AUCFeatRank(Y,L, AUCparam)
global VERBOSE

[~,nY] = size(Y);
D = zeros(nY,1);
if length(unique(L))==2
    uL = unique(L); 
    Ltemp = L; 
    Ltemp(L == min(uL)) = -1; 
    Ltemp(L ~= min(uL)) = 1; 
else 
    Ltemp = L; 
end
for i=1:nY
  D(i) = AUC(Ltemp,Y(:,i));  
end
if exist('AUCparam','var')
    ind = D<AUCparam;
    D=D-AUCparam;
    D(ind)=0;
    if VERBOSE, fprintf('\nAUC: Removed %g features below delta of %g', sum(ind),AUCparam); end
end    
if VERBOSE, fprintf('\nAUC: min = %g; max = %g', min(D), max(D)); end
end