function Y_mapped = nkdev_MapModelPredictionsMedian(O, R, I)

[m, n] = size(I);
Y_mapped = zeros(1,n);
%mR = nm_nanmean(R);
%sR = nm_nanstd(R);
for i=1:m
    Y_mapped(I(i,:)) = Y_mapped(I(i,:)) + R(i)-O;
end
Y_mapped = Y_mapped ./ (m * (sum(I)/m));
Y_mapped(isnan(Y_mapped)) = 0;

