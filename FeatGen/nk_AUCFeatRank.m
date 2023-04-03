function D = nk_AUCFeatRank(Y,L, AUCparam)
global VERBOSE

[~,nY] = size(Y);
D = zeros(nY,1);
for i=1:nY
  D(i) = AUC(L,Y(:,i));  
end
if exist('AUCparam','var')
    ind = D<AUCparam;
    D=D-AUCparam;
    D(ind)=0;
    if VERBOSE, fprintf('\nAUC: Removed %g features below delta of %g', sum(ind),AUCparam); end
end    
if VERBOSE, fprintf('\nAUC: min = %g; max = %g', min(D), max(D)); end
end