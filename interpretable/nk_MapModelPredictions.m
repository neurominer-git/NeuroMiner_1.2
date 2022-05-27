function Y_mapped = nk_MapModelPredictions(O, R, I, method)

[m, n] = size(I);
Y_mapped = zeros(1,n);
switch method
    case 'posneg'
        RangeD = abs(R(:,1) - R(:,2));
    case 'median'
        RangeD = R-O;
end
for i=1:m
    Y_mapped(I(i,:)) = Y_mapped(I(i,:)) + RangeD(i);
end
Y_mapped = Y_mapped./ sum(I);
Y_mapped(isnan(Y_mapped)) = 0;

