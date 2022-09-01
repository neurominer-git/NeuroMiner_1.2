function f = cal_mi(I1,I2)
% A function with the simplest form to
% compute Mutual Information of two images
% written by HU Xiubing
% Wuhan University,China
%  huxb@whu.edu.cn
size1 = size(I1);
size2 = size(I2);
if  ne(size1,size2)
    error('Please enter the same size images ! ');
end
I1 = I1(:);
I2 = I2(:);
min1 = min(I1); max1 = max(I1);
min2 = min(I2); max2 = max(I2);
N1 = max1 - min1 + 1;  % Grayscale of the image I1
N2 = max2 - min2 + 1;  % Grayscale of the image I2

ht = zeros(N1,N2);     % Joint histogram
for n = 1:length(I1);
    ht(I1(n)-min1+1,I2(n)-min2+1) = ht(I1(n)-min1+1,I2(n)-min2+1) + 1;
end

ht = ht/length(I1);  % normalized joint histogram
ym = sum(ht );       % sum of the rows of normalized joint histogram
xm = sum(ht');       % sum of columns of normalized joint histogran

Hy   =      sum(ym.*log2(ym+(ym==0)));
Hx   =      sum(xm.*log2(xm+(xm==0)));
h_xy =  sum(sum(ht.*log2(ht+(ht==0)))); % joint entropy

f = -(Hx+Hy)/h_xy;   % Mutual information
end