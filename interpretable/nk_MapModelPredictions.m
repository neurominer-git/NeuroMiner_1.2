function [Y_mapped, Y_mapped_ciu, Y_mapped_cil, Y_mapped_std] = nk_MapModelPredictions(n, O, R, I, MapIdx, Method, NormRange, ZnormData)

m = size(R,1);
if ~exist('ZnormData','var') || isempty('ZnormData')
    ZnormData = false;
end
tY_mapped = zeros(m,n);
switch Method
    case 'posneg'
        RangeD = (R(:,1) - R(:,2))-O;
    case 'median'
        RangeD = R-O;
end
for i=1:m
    tY_mapped(i,MapIdx(I(i,:))) = RangeD(i)*100/NormRange;
end
Y_mapped = zeros(1,n);
Y_mapped_ciu = zeros(1,n);
Y_mapped_cil = zeros(1,n);
Y_mapped_std = zeros(1,n);
for i=1:n
    idx = isfinite(tY_mapped(:,i)) & tY_mapped(:,i)~=0;
    ci = nm_95confint(tY_mapped(idx,i));
    Y_mapped_cil(i) = ci(1,:);
    Y_mapped_ciu(i) = ci(2,:);
    Y_mapped_std(i) = nm_nanstd(tY_mapped(idx,i));
    Y_mapped(i) = nm_nanmean(tY_mapped(idx,i));
end
switch ZnormData 
    case 2
        % mean centering 
        idx = Y_mapped ~=0 | isfinite(Y_mapped);
        mYmapped = nm_nanmean(Y_mapped(idx));
        Y_mapped(idx) = Y_mapped(idx) -  mYmapped ; 
    case 3
        % z-normalisation
        idx = Y_mapped ~=0 | isfinite(Y_mapped);
        mYmapped = nm_nanmean(Y_mapped(idx));
        sdYmapped = nm_nanmean(Y_mapped(idx));
        Y_mapped(idx) = ( Y_mapped(idx) -  mYmapped ) / sdYmapped;
        Y_mapped_cil(idx) = ( Y_mapped_cil(idx) -  mYmapped ) / sdYmapped;
        Y_mapped_ciu(idx) = ( Y_mapped_ciu(idx) -  mYmapped ) / sdYmapped;
        Y_mapped_std(idx) = ( Y_mapped_std(idx) -  mYmapped ) / sdYmapped;
end
idx = ~isfinite(Y_mapped);
Y_mapped(idx) = 0;
Y_mapped_cil(idx)=0;
Y_mapped_ciu(idx)=0;
Y_mapped_std(idx)=0;
