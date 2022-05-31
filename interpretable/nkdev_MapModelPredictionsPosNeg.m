function Y_mapped = nkdev_MapModelPredictionsPosNeg(O, R, I)

[m, n] = size(I);
Y_mapped = zeros(1,n);
RangeD = abs(R(:,1) - R(:,2));
for i=1:m
    Y_mapped(I(i,:)) = Y_mapped(I(i,:)) + RangeD(i);
end
Y_mapped = Y_mapped./ sum(I);
Y_mapped(isnan(Y_mapped)) = 0;

