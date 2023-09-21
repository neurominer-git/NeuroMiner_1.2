function [Y_mapped, Y_mapped_ciu, Y_mapped_cil, Y_mapped_std] = nk_MapModelPredictions(n, O, R, I, MapIdx, Method, NormRange, ZnormData)

m = size(R,1);
if ~exist('ZnormData','var') || isempty('ZnormData')
    ZnormData = false;
end
tY_mapped       = zeros(m,n);
Y_mapped        = zeros(1,n);
Y_mapped_ciu    = zeros(1,n);
Y_mapped_cil    = zeros(1,n);
Y_mapped_std    = zeros(1,n);

switch Method
    case 'posneg'
        RangeD = (R(:,1) - R(:,2))-O;
    case {'median','medianflip','random'}
        RangeD = R-O;
end
for i=1:m
    tY_mapped(i,MapIdx(I(i,:))) = RangeD(i)*100/NormRange;
end

switch ZnormData 
    case 2
        % mean centering 
        idx = tY_mapped ~=0 & isfinite(tY_mapped);
        mYmapped = nm_nanmean(tY_mapped(idx));
        tY_mapped(idx) = tY_mapped(idx) -  mYmapped ; 
    case 3
        % z-normalisation
        idx = tY_mapped ~=0 & isfinite(tY_mapped);
        mYmapped = nm_nanmean(tY_mapped(idx));
        sdYmapped = nm_nanstd(tY_mapped(idx));
        tY_mapped(idx) = ( tY_mapped(idx) -  mYmapped ) / sdYmapped;
end

for i=1:n
    idx = tY_mapped(:,i)~=0 & isfinite(tY_mapped(:,i)) ;
    if isempty(idx), continue; end
    ci = nm_95confint(tY_mapped(idx,i));
    Y_mapped_cil(i) = ci(1,:);
    Y_mapped_ciu(i) = ci(2,:);
    Y_mapped_std(i) = std(tY_mapped(idx,i));
    Y_mapped(i) = median(tY_mapped(idx,i));
end

