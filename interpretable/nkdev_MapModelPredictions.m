function Y_mapped = nkdev_MapModelPredictions(O, R, I)

[m, n] = size(I);
Y_mapped = zeros(1,n);
RangeD = R(:,1) - R(:,2);
RangeD = (RangeD-nm_nanmean(RangeD))*100/abs(O);
for i=1:m
    Y_mapped(I(i,:)) = Y_mapped(I(i,:)) + RangeD(i);
end
Y_mapped = Y_mapped./ sum(I);
Y_mapped(isnan(Y_mapped)) = 0;

